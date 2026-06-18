use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use crate::cluster::SeqDb;
use crate::cluster::kmer::extract_seeds;
use crate::error::{AppError, ErrorKind, Result};
use crate::index::bloom::BloomFilter;
use crate::index::codec::simd_dispatch::preferred_codec_for_cold;
use crate::index::codec::{CodecKind, delta_encode};
use crate::index::learned::PgmIndexV2;
use crate::index::manifest::{Manifest, ManifestSegment, next_segment_id, save_manifest};
use crate::index::segment::{SegmentMeta, segment_dir, write_codec_meta, write_segment_meta};
use crate::index::tier::{HotMode, TierAssignment, assign_tiers};
use crate::index::tune::pruning::{SegmentStats, write_segment_stats};
use crate::index::tune::seed_policy::{SeedPolicy, SeedPolicyMode, write_seed_policy};
use crate::io::atomic::write_atomic;

#[derive(Debug, Clone)]
pub struct BuildIndexConfig {
    pub dbdir: PathBuf,
    pub indexdir: PathBuf,
    pub k: usize,
    pub m: usize,
    pub seed: u64,
    pub hot_df_threshold: Option<u32>,
    pub hot_mode: HotMode,
    pub legacy_uncompressed: bool,
    pub pgm_epsilon: usize,
    pub seed_policy_mode: Option<SeedPolicyMode>,
}

pub fn build_fresh_index(cfg: &BuildIndexConfig) -> Result<u32> {
    fs::create_dir_all(&cfg.indexdir)
        .map_err(|e| AppError::io(format!("create dir {}", cfg.indexdir.display()), e))?;
    let id = 0;
    build_segment_with_global_base(cfg, id, 0)?;
    write_seed_policy(
        &cfg.indexdir,
        &SeedPolicy {
            mode: cfg
                .seed_policy_mode
                .clone()
                .unwrap_or(SeedPolicyMode::FixedK),
            k_values: vec![
                (cfg.k.saturating_sub(1)) as u32,
                cfg.k as u32,
                (cfg.k + 1) as u32,
            ],
            m: cfg.m as u32,
            rules: crate::index::tune::seed_policy::SeedPolicyRules {
                short_query_use_k: cfg.k.saturating_sub(1) as u32,
                long_query_use_k: cfg.k as u32,
            },
        },
    )?;
    let manifest = Manifest {
        version: 1,
        index_name: "kira-index".to_string(),
        segments: vec![ManifestSegment { segment_id: id }],
    };
    save_manifest(&cfg.indexdir, &manifest)?;
    Ok(id)
}

pub fn build_next_segment(
    cfg: &BuildIndexConfig,
    manifest: &Manifest,
    global_seq_base: u64,
) -> Result<u32> {
    let id = next_segment_id(manifest);
    build_segment_with_global_base(cfg, id, global_seq_base)?;
    Ok(id)
}

pub fn build_segment(cfg: &BuildIndexConfig, segment_id: u32) -> Result<()> {
    build_segment_with_global_base(cfg, segment_id, 0)
}

fn build_segment_with_global_base(
    cfg: &BuildIndexConfig,
    segment_id: u32,
    global_seq_base: u64,
) -> Result<()> {
    if cfg.k == 0 || cfg.m == 0 {
        return Err(AppError::new(ErrorKind::Validation, "k and m must be >= 1"));
    }

    if cfg.hot_mode == HotMode::Roaring && !cfg!(feature = "roaring") {
        return Err(AppError::new(
            ErrorKind::Unsupported,
            "ROARING_FEATURE_NOT_ENABLED",
        ));
    }

    let db = SeqDb::load(&cfg.dbdir)?;
    let mut pairs: Vec<(u64, u32)> = Vec::new();

    for seq_id in 0..db.n_seqs() as u32 {
        let seeds = extract_seeds(db.seq(seq_id), seq_id, cfg.k, cfg.m);
        for s in seeds {
            pairs.push((s.hash ^ cfg.seed, s.seq_id));
        }
    }

    pairs.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    let mut postings_by_key: BTreeMap<u64, Vec<u32>> = BTreeMap::new();
    for (key, seq_id) in pairs {
        let entry = postings_by_key.entry(key).or_default();
        if entry.last().copied() != Some(seq_id) {
            entry.push(seq_id);
        }
    }

    let mut keys = Vec::with_capacity(postings_by_key.len());
    let mut postings_lists = Vec::with_capacity(postings_by_key.len());
    let mut dfs = Vec::with_capacity(postings_by_key.len());
    for (key, list) in postings_by_key {
        dfs.push(list.len() as u32);
        keys.push(key);
        postings_lists.push(list);
    }

    let hot_df_threshold = cfg
        .hot_df_threshold
        .unwrap_or(((db.n_seqs() as u32) / 64).max(32));
    let tiers = assign_tiers(&dfs, hot_df_threshold);

    let seg_dir = segment_dir(&cfg.indexdir, segment_id);
    fs::create_dir_all(&seg_dir)
        .map_err(|e| AppError::io(format!("create dir {}", seg_dir.display()), e))?;

    write_atomic(&seg_dir.join("keys.u64"), &u64_to_bytes(&keys))?;

    let pgm = PgmIndexV2::build(&keys, cfg.pgm_epsilon.max(1))?;
    pgm.save_to_segment(&seg_dir)?;

    if cfg.legacy_uncompressed {
        write_legacy_postings(&seg_dir, &postings_lists)?;
    } else {
        write_compressed_postings(&seg_dir, &postings_lists, &dfs, &tiers, cfg.hot_mode)?;
        write_atomic(&seg_dir.join("df.u32"), &u32_to_bytes(&dfs))?;

        let bloom = BloomFilter::build(&keys, 10, 7, cfg.seed ^ 0xB100B100);
        write_atomic(&seg_dir.join("bloom.bin"), bloom.bits())?;
        let bloom_meta = serde_json::to_vec_pretty(bloom.meta()).map_err(|e| {
            AppError::new(ErrorKind::Internal, format!("serialize bloom meta: {e}"))
        })?;
        write_atomic(&seg_dir.join("bloom_meta.json"), &bloom_meta)?;
    }

    let stats = SegmentStats::from_data(segment_id, db.n_seqs() as u64, &keys, &dfs);
    write_segment_stats(&seg_dir, &stats)?;

    let db_fp = db_fingerprint(&cfg.dbdir)?;
    let meta = SegmentMeta {
        version: 2,
        segment_id,
        db_fingerprint: db_fp,
        db_path: cfg.dbdir.to_string_lossy().to_string(),
        k: cfg.k,
        m: cfg.m,
        seed: cfg.seed,
        n_seqs: db.n_seqs() as u64,
        global_seq_base,
        compressed: !cfg.legacy_uncompressed,
        hot_df_threshold,
        hot_mode: cfg.hot_mode.as_str().to_string(),
        codec_hot: if cfg.hot_mode == HotMode::Skip {
            "skip".to_string()
        } else {
            CodecKind::Bp128.name().to_string()
        },
        codec_cold: preferred_codec_for_cold().name().to_string(),
    };
    write_segment_meta(&seg_dir, &meta)?;

    Ok(())
}

fn write_legacy_postings(seg_dir: &Path, postings_lists: &[Vec<u32>]) -> Result<()> {
    let mut offsets = Vec::with_capacity(postings_lists.len() + 1);
    let mut postings = Vec::new();
    offsets.push(0u64);
    for list in postings_lists {
        postings.extend_from_slice(list);
        offsets.push(postings.len() as u64);
    }
    write_atomic(&seg_dir.join("offsets.u64"), &u64_to_bytes(&offsets))?;
    write_atomic(&seg_dir.join("postings.u32"), &u32_to_bytes(&postings))
}

fn write_compressed_postings(
    seg_dir: &Path,
    postings_lists: &[Vec<u32>],
    _dfs: &[u32],
    tiers: &[TierAssignment],
    hot_mode: HotMode,
) -> Result<()> {
    let codec_hot = CodecKind::Bp128;
    let codec_cold = preferred_codec_for_cold();

    let mut tiers_u8 = Vec::with_capacity(tiers.len());
    let mut tier_pos = Vec::with_capacity(tiers.len());

    let mut offs_hot = vec![0u64];
    let mut offs_cold = vec![0u64];
    let mut post_hot = Vec::new();
    let mut post_cold = Vec::new();
    let mut hot_idx = 0u32;
    let mut cold_idx = 0u32;

    let mut gaps = Vec::new();
    for (i, list) in postings_lists.iter().enumerate() {
        delta_encode(list, &mut gaps);

        match tiers[i] {
            TierAssignment::Hot => {
                tiers_u8.push(1);
                tier_pos.push(hot_idx);
                hot_idx += 1;
                if hot_mode == HotMode::Skip {
                    offs_hot.push(post_hot.len() as u64);
                } else {
                    let start = post_hot.len();
                    codec_hot.encode(&gaps, &mut post_hot);
                    let _ = start;
                    offs_hot.push(post_hot.len() as u64);
                }
            }
            TierAssignment::Cold => {
                tiers_u8.push(0);
                tier_pos.push(cold_idx);
                cold_idx += 1;
                codec_cold.encode(&gaps, &mut post_cold);
                offs_cold.push(post_cold.len() as u64);
            }
        }
    }

    write_atomic(&seg_dir.join("tiers.u8"), &tiers_u8)?;
    write_atomic(&seg_dir.join("tier_pos.u32"), &u32_to_bytes(&tier_pos))?;
    write_atomic(&seg_dir.join("offsets_hot.u64"), &u64_to_bytes(&offs_hot))?;
    write_atomic(&seg_dir.join("offsets_cold.u64"), &u64_to_bytes(&offs_cold))?;
    write_atomic(&seg_dir.join("postings_hot.bin"), &post_hot)?;
    write_atomic(&seg_dir.join("postings_cold.bin"), &post_cold)?;

    write_codec_meta(seg_dir, "codec_hot.json", codec_hot)?;
    write_codec_meta(seg_dir, "codec_cold.json", codec_cold)?;

    Ok(())
}

fn db_fingerprint(dbdir: &Path) -> Result<u64> {
    let bytes = std::fs::read(dbdir.join("meta.json"))
        .map_err(|e| AppError::io(format!("read {}/meta.json", dbdir.display()), e))?;
    Ok(crate::simd::hash_bytes32(&bytes))
}

fn u64_to_bytes(v: &[u64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(v.len() * 8);
    for x in v {
        out.extend_from_slice(&x.to_le_bytes());
    }
    out
}

fn u32_to_bytes(v: &[u32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(v.len() * 4);
    for x in v {
        out.extend_from_slice(&x.to_le_bytes());
    }
    out
}
