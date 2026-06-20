use rayon::prelude::*;

use crate::alignment::AlignmentResult;
use crate::alignment::banded::smith_waterman_gotoh_banded_tb;
use crate::alignment::farrar::sw_score;
use crate::cascade::CascadeThresholds;
use crate::cluster::filter::{
    coverage_prefilter, coverage_satisfied, hamming_filter, ungapped_filter,
};
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

/// Candidate lists at least this long use the rayon pool for the prefilter
/// stage. Smaller lists run sequentially to avoid per-task overhead (most
/// representatives have only a handful of candidates).
const PAR_PREFILTER_THRESHOLD: usize = 512;

/// Maximum Smith-Waterman DP matrix size (|a| * |b|) aligned exactly. Pairs
/// above this (very long sequences) are skipped rather than merged, bounding
/// worst-case per-pair work. 16M cells covers proteins up to ~4000 residues.
const MAX_DP_CELLS: u64 = 16_000_000;

/// Fraction of `min_identity * min_len` used as the Farrar score pre-screen
/// floor. Deliberately conservative so the screen never rejects a pair the full
/// alignment would accept (verified against the no-screen result on UniRef50).
const SCORE_FLOOR_FACTOR: f32 = 1.0;

#[derive(Clone, Copy)]
enum Screen {
    HammingReject,
    UngappedReject,
    Pass,
    NeedsGapped,
}

pub fn evaluate_candidates(
    db: &SeqDb,
    cfg: &ClusterConfig,
    thresholds: CascadeThresholds,
    rep_id: u32,
    candidates: &[u32],
    pool: Option<&rayon::ThreadPool>,
) -> (Vec<PairDecision>, EvalStats) {
    let rep = db.seq(rep_id);
    let mut stats = EvalStats {
        candidate_pairs_evaluated: candidates.len() as u64,
        ..EvalStats::default()
    };

    // Stage 1: prefilter every candidate (hamming + coverage + ungapped
    // identity). The candidate (query) is `seq`; the representative (target) is
    // `rep`. Coverage is handled here via cov_mode, so ungapped_filter is asked
    // for identity only (min_coverage = 0).
    let screen = |&sid: &u32| -> Screen {
        let seq = db.seq(sid);
        let qlen = seq.len();
        let tlen = rep.len();
        let overlap = qlen.min(tlen) as f32;
        let allowed = ((1.0 - thresholds.prefilter_identity) * overlap)
            .floor()
            .max(0.0) as u32;
        if !hamming_filter(rep, seq, allowed, cfg.backend) {
            return Screen::HammingReject;
        }
        // Gate on the FINAL coverage, not the loosened prefilter coverage:
        // alignment can never raise coverage above the length-ratio ceiling
        // (aligned_len <= min(q, t)), so a pair that cannot reach final_coverage
        // by length alone can never pass the final gate. Pruning it here avoids
        // a wasted (and potentially huge) gapped alignment with zero recall loss.
        if !coverage_prefilter(qlen, tlen, thresholds.final_coverage, cfg.cov_mode) {
            return Screen::UngappedReject;
        }
        if !ungapped_filter(rep, seq, thresholds.prefilter_identity, 0.0, cfg.backend) {
            return Screen::UngappedReject;
        }
        if thresholds.run_gapped {
            Screen::NeedsGapped
        } else {
            Screen::Pass
        }
    };

    let prefilter_t0 = std::time::Instant::now();
    let screens: Vec<Screen> = match pool {
        Some(p) if candidates.len() >= PAR_PREFILTER_THRESHOLD => {
            p.install(|| candidates.par_iter().map(screen).collect())
        }
        _ => candidates.iter().map(screen).collect(),
    };
    if let Some(p) = &cfg.profiler {
        p.add_stage_time("filter_prefilter", prefilter_t0.elapsed().as_nanos());
    }

    let mut decisions: Vec<PairDecision> = candidates
        .iter()
        .map(|&sid| PairDecision {
            seq_id: sid,
            passed: false,
        })
        .collect();
    let mut gapped_pairs: Vec<(u32, u32)> = Vec::new();
    let mut gapped_index: Vec<usize> = Vec::new();
    for (i, screen) in screens.iter().enumerate() {
        match screen {
            Screen::HammingReject => stats.hamming_rejected += 1,
            Screen::UngappedReject => stats.ungapped_rejected += 1,
            Screen::Pass => decisions[i].passed = true,
            Screen::NeedsGapped => {
                gapped_pairs.push((rep_id, candidates[i]));
                gapped_index.push(i);
            }
        }
    }

    if gapped_pairs.is_empty() {
        return (decisions, stats);
    }

    let gapped_t0 = std::time::Instant::now();
    let gapped_results = run_gapped_stage(db, cfg, &gapped_pairs, pool);
    if let Some(p) = &cfg.profiler {
        p.add_stage_time("filter_gapped", gapped_t0.elapsed().as_nanos());
    }

    for (local_idx, aln) in gapped_results.into_iter().enumerate() {
        let decision_idx = gapped_index[local_idx];
        if aln.aligned_len == 0 {
            stats.gapped_rejected += 1;
            continue;
        }
        let seq = db.seq(candidates[decision_idx]);
        let identity = (aln.matches as f32) / (aln.aligned_len as f32);
        let coverage_ok = coverage_satisfied(
            aln.aligned_len,
            seq.len(),
            rep.len(),
            thresholds.final_coverage,
            cfg.cov_mode,
        );
        let passed = identity >= thresholds.final_identity && coverage_ok;
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

    let align_pair = |a: u32, b: u32| -> AlignmentResult {
        let sa = db.seq(a);
        let sb = db.seq(b);
        // Skip pathologically large pairs (this data set has a 49k-residue
        // protein): even a banded alignment over such a pair is wasteful, and a
        // pair with such mismatched length cannot pass the coverage gate anyway.
        // Counted as a gapped rejection (not merged).
        if (sa.len() as u64).saturating_mul(sb.len() as u64) > MAX_DP_CELLS {
            return AlignmentResult {
                score: 0,
                aligned_len: 0,
                matches: 0,
            };
        }
        // Fast vectorized score pre-screen (Farrar): a pair whose full optimal
        // score is far below what an identity>=min_id, coverage>=cov alignment
        // needs cannot pass, so skip the (more expensive) traceback DP. The full
        // SW score is an upper bound on the banded score, so this never drops a
        // pair the banded DP would have accepted on score grounds.
        let min_len = sa.len().min(sb.len());
        let score_floor = (min_len as f32 * cfg.min_identity * SCORE_FLOOR_FACTOR) as i32;
        if sw_score(sa, sb, &cfg.scoring, cfg.backend) < score_floor {
            return AlignmentResult {
                score: 0,
                aligned_len: 0,
                matches: 0,
            };
        }

        // Banded row-DP with traceback for exact matches/aligned_len.
        let half_band = band_half_width(sa.len(), sb.len());
        smith_waterman_gotoh_banded_tb(sa, sb, &cfg.scoring, half_band)
    };

    match pool {
        Some(p) => p.install(|| pairs.par_iter().map(|&(a, b)| align_pair(a, b)).collect()),
        None => pairs.iter().map(|&(a, b)| align_pair(a, b)).collect(),
    }
}

/// Half-width of the alignment band around the start→end diagonal. The band
/// absorbs internal indel wiggle (the length difference is already covered by
/// the diagonal slope). 1/4 of the longer sequence reproduced the full
/// (unbanded) clustering result bit-for-bit on UniRef50; 1/8 was ~1.6x faster
/// but lost a handful of merges, so 1/4 is kept as the lossless default. Short
/// sequences get a band >= their length, i.e. a full alignment.
fn band_half_width(la: usize, lb: usize) -> usize {
    (la.max(lb) / 4).max(64)
}

pub fn gpu_batches_for_pairs(pairs: &[(u32, u32)], cfg: &ClusterConfig) -> Vec<batch::GpuBatch> {
    build_batches(pairs, &cfg.scoring, cfg.batch_size)
}
