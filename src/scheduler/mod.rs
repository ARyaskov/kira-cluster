use rayon::prelude::*;

use crate::alignment::{AlignmentResult, smith_waterman_gotoh};
use crate::cascade::CascadeThresholds;
use crate::cluster::filter::{hamming_filter, ungapped_filter};
use crate::cluster::{ClusterConfig, SeqDb};
use crate::gpu::{GpuConfig, align_pairs_on_gpu};
use crate::scheduler::batch::build_batches;

pub mod batch;

#[derive(Debug, Clone, Copy)]
pub struct PairDecision {
    pub seq_id: u32,
    pub passed: bool,
}

#[derive(Debug, Clone, Default)]
pub struct EvalStats {
    pub candidate_pairs_evaluated: u64,
    pub hamming_rejected: u64,
    pub ungapped_rejected: u64,
    pub gapped_rejected: u64,
}

pub fn evaluate_candidates(
    db: &SeqDb,
    cfg: &ClusterConfig,
    thresholds: CascadeThresholds,
    rep_id: u32,
    candidates: &[u32],
    pool: Option<&rayon::ThreadPool>,
) -> (Vec<PairDecision>, EvalStats) {
    let mut decisions: Vec<PairDecision> = candidates
        .iter()
        .map(|&sid| PairDecision {
            seq_id: sid,
            passed: false,
        })
        .collect();

    let mut gapped_pairs: Vec<(u32, u32)> = Vec::new();
    let mut gapped_index: Vec<usize> = Vec::new();

    let mut stats = EvalStats {
        candidate_pairs_evaluated: candidates.len() as u64,
        ..EvalStats::default()
    };
    let mut hamming_ns = 0u128;
    let mut ungapped_ns = 0u128;

    for (i, &sid) in candidates.iter().enumerate() {
        let rep = db.seq(rep_id);
        let seq = db.seq(sid);
        let overlap = rep.len().min(seq.len()) as f32;
        let allowed = ((1.0 - thresholds.prefilter_identity) * overlap)
            .floor()
            .max(0.0) as u32;

        let t0 = std::time::Instant::now();
        let ham_ok = hamming_filter(rep, seq, allowed, cfg.backend);
        hamming_ns += t0.elapsed().as_nanos();
        if !ham_ok {
            stats.hamming_rejected += 1;
            continue;
        }
        let t1 = std::time::Instant::now();
        let ungapped_ok = ungapped_filter(
            rep,
            seq,
            thresholds.prefilter_identity,
            thresholds.prefilter_coverage,
            cfg.backend,
        );
        ungapped_ns += t1.elapsed().as_nanos();
        if !ungapped_ok {
            stats.ungapped_rejected += 1;
            continue;
        }

        if !thresholds.run_gapped {
            decisions[i].passed = true;
            continue;
        }

        gapped_pairs.push((rep_id, sid));
        gapped_index.push(i);
    }

    if let Some(p) = &cfg.profiler {
        p.add_stage_time("filter_hamming", hamming_ns);
        p.add_stage_time("filter_ungapped", ungapped_ns);
    }

    if gapped_pairs.is_empty() {
        return (decisions, stats);
    }

    let gapped_t0 = std::time::Instant::now();
    let gapped_results = run_gapped_stage(db, cfg, thresholds, &gapped_pairs, pool);
    let gapped_ns = gapped_t0.elapsed().as_nanos();
    if let Some(p) = &cfg.profiler {
        p.add_stage_time("filter_gapped", gapped_ns);
    }

    for (local_idx, aln) in gapped_results.into_iter().enumerate() {
        let decision_idx = gapped_index[local_idx];
        if aln.aligned_len == 0 {
            stats.gapped_rejected += 1;
            continue;
        }
        let rep = db.seq(rep_id);
        let seq = db.seq(candidates[decision_idx]);
        let identity = (aln.matches as f32) / (aln.aligned_len as f32);
        let coverage = (aln.aligned_len as f32) / (rep.len().min(seq.len()) as f32);

        let passed = identity >= thresholds.final_identity && coverage >= thresholds.final_coverage;
        decisions[decision_idx].passed = passed;
        if !passed {
            stats.gapped_rejected += 1;
        }
    }

    (decisions, stats)
}

fn run_gapped_stage(
    db: &SeqDb,
    cfg: &ClusterConfig,
    thresholds: CascadeThresholds,
    pairs: &[(u32, u32)],
    pool: Option<&rayon::ThreadPool>,
) -> Vec<AlignmentResult> {
    if cfg.use_gpu {
        let gpu_cfg = GpuConfig {
            backend: cfg.gpu_backend,
            batch_size: cfg.batch_size,
            gpu_memory_limit: cfg.gpu_memory_limit,
        };
        if let Ok(results) = align_pairs_on_gpu(db, pairs, &cfg.scoring, gpu_cfg) {
            return results;
        }
    }

    let run_cpu = || {
        pairs
            .iter()
            .map(|&(a, b)| {
                let sa = db.seq(a);
                let sb = db.seq(b);
                if thresholds.short_gapped {
                    let (xa, xb) = short_views(sa, sb);
                    smith_waterman_gotoh(xa, xb, &cfg.scoring, cfg.backend)
                } else {
                    smith_waterman_gotoh(sa, sb, &cfg.scoring, cfg.backend)
                }
            })
            .collect::<Vec<_>>()
    };

    match pool {
        Some(p) => p.install(|| {
            pairs
                .par_iter()
                .map(|&(a, b)| {
                    let sa = db.seq(a);
                    let sb = db.seq(b);
                    if thresholds.short_gapped {
                        let (xa, xb) = short_views(sa, sb);
                        smith_waterman_gotoh(xa, xb, &cfg.scoring, cfg.backend)
                    } else {
                        smith_waterman_gotoh(sa, sb, &cfg.scoring, cfg.backend)
                    }
                })
                .collect::<Vec<_>>()
        }),
        None => run_cpu(),
    }
}

fn short_views<'a>(a: &'a [u8], b: &'a [u8]) -> (&'a [u8], &'a [u8]) {
    const BAND_LEN: usize = 256;
    (&a[..a.len().min(BAND_LEN)], &b[..b.len().min(BAND_LEN)])
}

pub fn gpu_batches_for_pairs(pairs: &[(u32, u32)], cfg: &ClusterConfig) -> Vec<batch::GpuBatch> {
    build_batches(pairs, &cfg.scoring, cfg.batch_size)
}
