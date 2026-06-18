use std::path::{Path, PathBuf};

use memmap2::Mmap;
use serde::{Deserialize, Serialize};

use crate::cluster::SeqDb;
use crate::error::{AppError, ErrorKind, Result};
use crate::index::bloom::{BloomFilter, BloomMeta};
use crate::index::codec::simd_dispatch::decode_bp128_simd;
use crate::index::codec::{CodecKind, delta_decode};
use crate::index::learned::PgmIndexV2;
use crate::index::tier::HotMode;
use crate::io::atomic::write_atomic;
use crate::io::mmap::map_readonly;
use crate::simd::SimdBackend;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SegmentMeta {
    pub version: u32,
    pub segment_id: u32,
    pub db_fingerprint: u64,
    pub db_path: String,
    pub k: usize,
    pub m: usize,
    pub seed: u64,
    pub n_seqs: u64,
    #[serde(default)]
    pub global_seq_base: u64,
    #[serde(default)]
    pub compressed: bool,
    #[serde(default)]
    pub hot_df_threshold: u32,
    #[serde(default)]
    pub hot_mode: String,
    #[serde(default)]
    pub codec_hot: String,
    #[serde(default)]
    pub codec_cold: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CodecMeta {
    pub codec: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SegmentSeqRef<'a> {
    pub global_id: u64,
    pub segment_id: u32,
    pub local_id: u32,
    pub name: &'a str,
}

enum SegmentStorage {
    Legacy {
        offsets_map: Mmap,
        postings_map: Mmap,
    },
    Tiered {
        df_map: Mmap,
        tiers_map: Mmap,
        tier_pos_map: Mmap,
        offsets_hot_map: Mmap,
        offsets_cold_map: Mmap,
        postings_hot_map: Mmap,
        postings_cold_map: Mmap,
        codec_hot: CodecKind,
        codec_cold: CodecKind,
        hot_mode: HotMode,
        bloom: Option<BloomFilter>,
    },
}

pub struct SegmentHandle {
    pub meta: SegmentMeta,
    pub db: SeqDb,
    keys_map: Mmap,
    pgm: PgmIndexV2,
    storage: SegmentStorage,
}

impl SegmentHandle {
    pub fn open(indexdir: &Path, segment_id: u32) -> Result<Self> {
        let dir = segment_dir(indexdir, segment_id);
        let meta_map = map_readonly(&dir.join("segment_meta.json"))?;
        let meta: SegmentMeta = serde_json::from_slice(&meta_map)
            .map_err(|e| AppError::new(ErrorKind::Parse, format!("parse segment meta: {e}")))?;

        let keys_map = map_readonly(&dir.join("keys.u64"))?;
        let keys_len = keys_map.len() / 8;
        let pgm = if dir.join("kv_pgm.bin").exists() {
            PgmIndexV2::load_from_segment(&dir)?
        } else {
            PgmIndexV2::build(as_u64_slice(&keys_map), 1024.max(keys_len / 64))?
        };
        let db = SeqDb::load(Path::new(&meta.db_path))?;

        let storage = if dir.join("offsets_hot.u64").exists() {
            let codec_hot = read_codec_kind(&dir.join("codec_hot.json"))?;
            let codec_cold = read_codec_kind(&dir.join("codec_cold.json"))?;
            let hot_mode = match meta.hot_mode.as_str() {
                "skip" => HotMode::Skip,
                "roaring" => HotMode::Roaring,
                _ => HotMode::Bp128,
            };
            let bloom = if dir.join("bloom_meta.json").exists() && dir.join("bloom.bin").exists() {
                let mm = map_readonly(&dir.join("bloom_meta.json"))?;
                let meta: BloomMeta = serde_json::from_slice(&mm).map_err(|e| {
                    AppError::new(ErrorKind::Parse, format!("parse bloom meta: {e}"))
                })?;
                let bits = std::fs::read(dir.join("bloom.bin"))
                    .map_err(|e| AppError::io(format!("read {}/bloom.bin", dir.display()), e))?;
                Some(BloomFilter::from_parts(meta, bits))
            } else {
                None
            };

            SegmentStorage::Tiered {
                df_map: map_readonly(&dir.join("df.u32"))?,
                tiers_map: map_readonly(&dir.join("tiers.u8"))?,
                tier_pos_map: map_readonly(&dir.join("tier_pos.u32"))?,
                offsets_hot_map: map_readonly(&dir.join("offsets_hot.u64"))?,
                offsets_cold_map: map_readonly(&dir.join("offsets_cold.u64"))?,
                postings_hot_map: map_readonly(&dir.join("postings_hot.bin"))?,
                postings_cold_map: map_readonly(&dir.join("postings_cold.bin"))?,
                codec_hot,
                codec_cold,
                hot_mode,
                bloom,
            }
        } else {
            SegmentStorage::Legacy {
                offsets_map: map_readonly(&dir.join("offsets.u64"))?,
                postings_map: map_readonly(&dir.join("postings.u32"))?,
            }
        };

        Ok(Self {
            meta,
            db,
            keys_map,
            pgm,
            storage,
        })
    }

    pub fn keys(&self) -> &[u64] {
        as_u64_slice(&self.keys_map)
    }

    pub fn codec_hot_name(&self) -> &'static str {
        match &self.storage {
            SegmentStorage::Legacy { .. } => "legacy_u32",
            SegmentStorage::Tiered { codec_hot, .. } => codec_hot.name(),
        }
    }

    pub fn codec_cold_name(&self) -> &'static str {
        match &self.storage {
            SegmentStorage::Legacy { .. } => "legacy_u32",
            SegmentStorage::Tiered { codec_cold, .. } => codec_cold.name(),
        }
    }

    pub fn find_key(&self, key: u64) -> Option<usize> {
        if let SegmentStorage::Tiered { bloom: Some(b), .. } = &self.storage {
            if !b.maybe_contains(key) {
                return None;
            }
        }

        let keys = self.keys();
        if let Some(idx) = self.pgm.find_key(key) {
            if keys.get(idx).copied() == Some(key) {
                return Some(idx);
            }
        }

        keys.binary_search(&key).ok()
    }

    pub fn global_seq_id(&self, local_seq_id: u32) -> u64 {
        self.meta.global_seq_base + local_seq_id as u64
    }

    pub fn local_seq_id(&self, global_seq_id: u64) -> Option<u32> {
        let rel = global_seq_id.checked_sub(self.meta.global_seq_base)?;
        if rel < self.meta.n_seqs {
            u32::try_from(rel).ok()
        } else {
            None
        }
    }

    pub fn seq_ref(&self, local_seq_id: u32) -> SegmentSeqRef<'_> {
        SegmentSeqRef {
            global_id: self.global_seq_id(local_seq_id),
            segment_id: self.meta.segment_id,
            local_id: local_seq_id,
            name: self.db.name(local_seq_id),
        }
    }

    pub fn resolve_global_seq_id(&self, global_seq_id: u64) -> Option<SegmentSeqRef<'_>> {
        self.local_seq_id(global_seq_id)
            .map(|local_id| self.seq_ref(local_id))
    }

    pub fn df(&self, idx: usize) -> u32 {
        match &self.storage {
            SegmentStorage::Legacy {
                offsets_map,
                postings_map: _,
            } => {
                let offs = as_u64_slice(offsets_map);
                (offs[idx + 1] - offs[idx]) as u32
            }
            SegmentStorage::Tiered { df_map, .. } => as_u32_slice(df_map)[idx],
        }
    }

    pub fn postings_for_key_index(&self, idx: usize, backend: SimdBackend) -> Option<Vec<u32>> {
        match &self.storage {
            SegmentStorage::Legacy {
                offsets_map,
                postings_map,
            } => {
                let offs = as_u64_slice(offsets_map);
                let posts = as_u32_slice(postings_map);
                let start = offs[idx] as usize;
                let end = offs[idx + 1] as usize;
                Some(posts[start..end].to_vec())
            }
            SegmentStorage::Tiered {
                tiers_map,
                tier_pos_map,
                offsets_hot_map,
                offsets_cold_map,
                postings_hot_map,
                postings_cold_map,
                codec_hot,
                codec_cold,
                hot_mode,
                ..
            } => {
                let tiers = tiers_map.as_ref();
                let tier_pos = as_u32_slice(tier_pos_map);
                let tier = tiers[idx];
                let tier_idx = tier_pos[idx] as usize;

                if tier == 1 {
                    if *hot_mode == HotMode::Skip {
                        return None;
                    }
                    let offs = as_u64_slice(offsets_hot_map);
                    let start = offs[tier_idx] as usize;
                    let end = offs[tier_idx + 1] as usize;
                    let encoded = &postings_hot_map[start..end];
                    Some(decode_postings(*codec_hot, encoded, backend))
                } else {
                    let offs = as_u64_slice(offsets_cold_map);
                    let start = offs[tier_idx] as usize;
                    let end = offs[tier_idx + 1] as usize;
                    let encoded = &postings_cold_map[start..end];
                    Some(decode_postings(*codec_cold, encoded, backend))
                }
            }
        }
    }
}

fn decode_postings(codec: CodecKind, encoded: &[u8], backend: SimdBackend) -> Vec<u32> {
    let mut gaps = Vec::new();
    if codec == CodecKind::Bp128 {
        decode_bp128_simd(encoded, &mut gaps, backend);
    } else {
        codec.decode(encoded, &mut gaps);
    }
    let mut ids = Vec::new();
    delta_decode(&gaps, &mut ids);
    ids
}

fn read_codec_kind(path: &Path) -> Result<CodecKind> {
    let mm = map_readonly(path)?;
    let meta: CodecMeta = serde_json::from_slice(&mm)
        .map_err(|e| AppError::new(ErrorKind::Parse, format!("parse {}: {e}", path.display())))?;
    CodecKind::from_name(&meta.codec)
        .ok_or_else(|| AppError::new(ErrorKind::Parse, format!("unknown codec {}", meta.codec)))
}

fn as_u64_slice(map: &Mmap) -> &[u64] {
    // SAFETY: backing bytes are generated by this program as little-endian u64 table.
    unsafe { std::slice::from_raw_parts(map.as_ptr() as *const u64, map.len() / 8) }
}

fn as_u32_slice(map: &Mmap) -> &[u32] {
    // SAFETY: backing bytes are generated by this program as little-endian u32 table.
    unsafe { std::slice::from_raw_parts(map.as_ptr() as *const u32, map.len() / 4) }
}

pub fn segment_dir(indexdir: &Path, segment_id: u32) -> PathBuf {
    indexdir.join(format!("segment_{segment_id:04}"))
}

pub fn write_segment_meta(dir: &Path, meta: &SegmentMeta) -> Result<()> {
    let bytes = serde_json::to_vec_pretty(meta)
        .map_err(|e| AppError::new(ErrorKind::Internal, format!("serialize segment meta: {e}")))?;
    write_atomic(&dir.join("segment_meta.json"), &bytes)
}

pub fn write_codec_meta(dir: &Path, name: &str, codec: CodecKind) -> Result<()> {
    let meta = CodecMeta {
        codec: codec.name().to_string(),
    };
    let bytes = serde_json::to_vec_pretty(&meta)
        .map_err(|e| AppError::new(ErrorKind::Internal, format!("serialize codec meta: {e}")))?;
    write_atomic(&dir.join(name), &bytes)
}
