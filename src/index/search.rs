use std::io::Write;
use std::path::Path;

use crate::alignment::scoring::ScoringParams;
use crate::cluster::kmer::{KmerSeed, extract_seeds};
use crate::error::{AppError, ErrorKind, Result};
use crate::index::IndexHandle;
use crate::index::df::idf;
use crate::index::scoring::score_candidate;
use crate::index::segment::SegmentHandle;
use crate::index::tune::seed_policy::{SeedPolicy, choose_k_for_query};
use crate::io::atomic::write_atomic_with;
use crate::seq::fasta;
use crate::simd::SimdBackend;

#[derive(Debug, Clone)]
pub struct SearchOpts {
    pub k: usize,
    pub m: usize,
    pub top_k: usize,
    pub min_seed_hits: u32,
    pub min_identity: f32,
    pub min_coverage: f32,
    pub max_df: Option<u32>,
    pub prune_df_quantile: f32,
    pub max_seeds_per_query: u32,
    pub work_budget: Option<u64>,
    pub scoring: ScoringParams,
    pub backend: SimdBackend,
    pub seed: u64,
    pub seed_policy: Option<SeedPolicy>,
}

#[derive(Debug, Clone)]
pub struct Hit {
    pub query_id: String,
    pub target_id: u64,
    pub target_name: String,
    pub local_id: u32,
    pub score: i32,
    pub identity: f32,
    pub coverage: f32,
    pub seeds_hit: u32,
    pub segment_id: u32,
    pub explain: String,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct SearchStats {
    pub queries: u64,
    pub queries_with_hits: u64,
    pub total_hits: u64,
}

impl SearchStats {
    pub fn zero_hit_queries(&self) -> u64 {
        self.queries.saturating_sub(self.queries_with_hits)
    }
}

/// Verify that the search-time k/seed match what the index was built with.
/// A mismatch otherwise produces a silently empty result file: the probe key is
/// `hash ^ seed` over k-mers of length k, so any divergence misses every key.
fn validate_search_opts(handle: &IndexHandle, opts: &SearchOpts) -> Result<()> {
    for seg in &handle.segments {
        if opts.k != seg.meta.k {
            return Err(AppError::new(
                ErrorKind::Validation,
                format!(
                    "search --k {} does not match index k {} (segment {}); rebuild the index or pass --k {}",
                    opts.k, seg.meta.k, seg.meta.segment_id, seg.meta.k
                ),
            ));
        }
        if opts.seed != seg.meta.seed {
            return Err(AppError::new(
                ErrorKind::Validation,
                format!(
                    "search --seed {} does not match index seed {} (segment {}); pass --seed {}",
                    opts.seed, seg.meta.seed, seg.meta.segment_id, seg.meta.seed
                ),
            ));
        }
        if opts.m != seg.meta.m {
            eprintln!(
                "WARN KW2001 SEARCH_M_MISMATCH query_m={} index_m={} segment={}",
                opts.m, seg.meta.m, seg.meta.segment_id
            );
        }
    }
    Ok(())
}

pub fn search_index(
    handle: &IndexHandle,
    query_fasta: &Path,
    out_tsv: &Path,
    opts: &SearchOpts,
) -> Result<SearchStats> {
    validate_search_opts(handle, opts)?;

    let mut stats = SearchStats::default();
    write_atomic_with(out_tsv, |w| {
        let mut line = String::new();
        fasta::parse_fasta(query_fasta, |query_name, query_seq| {
            let mut hits = Vec::new();
            search_query_sequence(&handle.segments, &query_seq, opts, &mut hits);

            stats.queries += 1;
            if !hits.is_empty() {
                stats.queries_with_hits += 1;
            }
            stats.total_hits += hits.len() as u64;

            for mut h in hits {
                h.query_id = query_name.clone();
                line.clear();
                use std::fmt::Write as _;
                let _ = writeln!(
                    line,
                    "{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}",
                    h.query_id,
                    h.target_id,
                    h.score,
                    h.identity,
                    h.coverage,
                    h.seeds_hit,
                    h.segment_id,
                    h.explain,
                    h.target_name,
                    h.local_id
                );
                w.write_all(line.as_bytes())
                    .map_err(|e| AppError::io("write hits tsv", e))?;
            }
            Ok(())
        })
    })?;

    Ok(stats)
}

pub fn search_query_sequence(
    segments: &[SegmentHandle],
    query: &[u8],
    opts: &SearchOpts,
    out_hits: &mut Vec<Hit>,
) {
    let effective_k = opts
        .seed_policy
        .as_ref()
        .map(|p| choose_k_for_query(p, opts.k, query))
        .unwrap_or(opts.k)
        .max(1);

    let seeds = extract_seeds(query, 0, effective_k, opts.m);
    let selected_seeds = adaptive_seed_subset(&seeds, segments, opts);

    for seg in segments {
        let mut raw_pairs: Vec<(u32, i32, u64)> = Vec::new();
        let mut seg_seed_hits = 0usize;
        let mut estimated_work = 0u64;

        for s in &selected_seeds {
            let key = s.hash ^ opts.seed;
            let Some(idx) = seg.find_key(key) else {
                continue;
            };

            seg_seed_hits += 1;
            let df = seg.df(idx);
            if let Some(max_df) = opts.max_df {
                if df > max_df {
                    continue;
                }
            }

            estimated_work += df as u64;
            let weight = (idf(seg.meta.n_seqs, df) * 1024.0).round() as i32;
            if let Some(postings) = seg.postings_for_key_index(idx, opts.backend) {
                for seq_id in postings {
                    raw_pairs.push((seq_id, weight.max(1), key));
                }
            }

            if let Some(budget) = opts.work_budget {
                if estimated_work >= budget && raw_pairs.len() as u64 >= budget {
                    break;
                }
            }
        }

        if seg_seed_hits == 0 {
            continue;
        }
        if let Some(budget) = opts.work_budget {
            if estimated_work > budget.saturating_mul(4)
                && seg_seed_hits <= (selected_seeds.len() / 10).max(1)
            {
                continue;
            }
        }

        raw_pairs.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.2.cmp(&b.2)).then(b.1.cmp(&a.1)));

        let mut candidates: Vec<(u32, u32, i32, Vec<u64>)> = Vec::new();
        let mut i = 0usize;
        while i < raw_pairs.len() {
            let seq_id = raw_pairs[i].0;
            let mut seed_hits = 0u32;
            let mut weight_sum = 0i32;
            let mut top_keys = Vec::new();
            while i < raw_pairs.len() && raw_pairs[i].0 == seq_id {
                seed_hits += 1;
                weight_sum += raw_pairs[i].1;
                if top_keys.len() < 4 && top_keys.last().copied() != Some(raw_pairs[i].2) {
                    top_keys.push(raw_pairs[i].2);
                }
                i += 1;
            }
            candidates.push((seq_id, seed_hits, weight_sum, top_keys));
        }

        let mut cand_hits = Vec::new();
        for (seq_id, seeds_hit, seed_weight, top_keys) in candidates {
            if seeds_hit < opts.min_seed_hits {
                continue;
            }
            if let Some(sc) = score_candidate(
                query,
                seg.db.seq(seq_id),
                seeds_hit,
                opts.min_identity,
                opts.min_coverage,
                &opts.scoring,
                opts.backend,
            ) {
                let final_score = sc.score + seed_weight / 256;
                let target = seg.seq_ref(seq_id);
                let explain = format!(
                    "k={};keys={:x?};count={};w={}",
                    effective_k, top_keys, seeds_hit, seed_weight
                );
                cand_hits.push(Hit {
                    query_id: String::new(),
                    target_id: target.global_id,
                    target_name: target.name.to_string(),
                    local_id: target.local_id,
                    score: final_score,
                    identity: sc.identity,
                    coverage: sc.coverage,
                    seeds_hit,
                    segment_id: seg.meta.segment_id,
                    explain,
                });
            }
        }

        cand_hits.sort_by(|a, b| {
            b.score
                .cmp(&a.score)
                .then(a.target_id.cmp(&b.target_id))
                .then(a.segment_id.cmp(&b.segment_id))
        });
        if cand_hits.len() > opts.top_k {
            cand_hits.truncate(opts.top_k);
        }
        out_hits.extend(cand_hits);
    }

    out_hits.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.target_id.cmp(&b.target_id))
            .then(a.segment_id.cmp(&b.segment_id))
    });
    if out_hits.len() > opts.top_k {
        out_hits.truncate(opts.top_k);
    }
}

fn adaptive_seed_subset(
    seeds: &[KmerSeed],
    segments: &[SegmentHandle],
    opts: &SearchOpts,
) -> Vec<KmerSeed> {
    if seeds.is_empty() {
        return Vec::new();
    }

    let mut seed_df: Vec<(KmerSeed, u32)> = Vec::with_capacity(seeds.len());
    for s in seeds {
        let mut best_df = u32::MAX;
        for seg in segments {
            let key = s.hash ^ opts.seed;
            if let Some(idx) = seg.find_key(key) {
                best_df = best_df.min(seg.df(idx));
            }
        }
        seed_df.push((*s, best_df));
    }

    let q = opts.prune_df_quantile.clamp(0.0, 1.0);
    let threshold = if seed_df.is_empty() {
        u32::MAX
    } else {
        let mut vals: Vec<u32> = seed_df.iter().map(|(_, df)| *df).collect();
        vals.sort_unstable();
        let idx = ((vals.len() - 1) as f32 * q).round() as usize;
        vals[idx.min(vals.len() - 1)]
    };

    let mut kept: Vec<(KmerSeed, u32)> = seed_df
        .into_iter()
        .filter(|(_, df)| *df <= threshold)
        .collect();

    if kept.is_empty() {
        kept = seeds.iter().copied().map(|s| (s, u32::MAX)).collect();
    }

    kept.sort_unstable_by(|a, b| {
        a.1.cmp(&b.1)
            .then(a.0.hash.cmp(&b.0.hash))
            .then(a.0.pos.cmp(&b.0.pos))
    });
    let cap = opts.max_seeds_per_query.max(1) as usize;
    if kept.len() > cap {
        kept.truncate(cap);
    }

    kept.into_iter().map(|(s, _)| s).collect()
}
