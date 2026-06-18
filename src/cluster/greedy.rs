use crate::cascade::thresholds_for;
use crate::cluster::table::KmerTable;
use crate::cluster::{ClusterConfig, SeqDb};
use crate::error::{AppError, ErrorKind, Result};
use crate::profile::WorkCounters;
use crate::scheduler::evaluate_candidates;

#[derive(Debug, Clone)]
pub struct ClusterResult {
    pub representatives: Vec<u32>,
    pub assignments: Vec<u32>,
}

pub fn cluster(db: &SeqDb, cfg: &ClusterConfig, k: usize) -> Result<ClusterResult> {
    let pool = if cfg.cpu_threads > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(cfg.cpu_threads)
                .build()
                .map_err(|e| {
                    AppError::new(ErrorKind::Internal, format!("build thread pool: {e}"))
                })?,
        )
    } else {
        None
    };
    let pool_ref = pool.as_ref();

    let table = KmerTable::build(db, k, cfg.kmer_per_seq, pool_ref, cfg.profiler.as_deref());
    let grouping_guard = cfg.profiler.as_ref().map(|p| p.stage("candidate_grouping"));
    let rep_candidates = build_group_candidates(db, &table);
    drop(grouping_guard);
    let thresholds = thresholds_for(
        cfg.cascade_level,
        cfg.min_identity,
        cfg.min_coverage,
        cfg.sensitivity,
    );

    let mut order: Vec<u32> = (0..db.n_seqs() as u32).collect();
    order.sort_by(|a, b| db.seq_len(*b).cmp(&db.seq_len(*a)).then(a.cmp(b)));

    let mut assignments: Vec<Option<u32>> = vec![None; db.n_seqs()];
    let mut representatives: Vec<u32> = Vec::new();

    let greedy_guard = cfg.profiler.as_ref().map(|p| p.stage("greedy_assignment"));
    let mut counters = WorkCounters::default();

    for rep_id in order {
        if assignments[rep_id as usize].is_some() {
            continue;
        }

        let cid = representatives.len() as u32;
        representatives.push(rep_id);
        assignments[rep_id as usize] = Some(cid);

        let posting_len = rep_candidates[rep_id as usize].len() as u64;
        let mut pending: Vec<u32> = rep_candidates[rep_id as usize]
            .iter()
            .copied()
            .filter(|sid| assignments[*sid as usize].is_none())
            .collect();
        pending.sort_unstable();
        pending.dedup();

        let (passes, eval_stats) =
            evaluate_candidates(db, cfg, thresholds, rep_id, &pending, pool_ref);
        counters.total_posting_length_read += posting_len;
        counters.candidate_pairs_evaluated += eval_stats.candidate_pairs_evaluated;
        counters.hamming_rejected += eval_stats.hamming_rejected;
        counters.ungapped_rejected += eval_stats.ungapped_rejected;
        counters.gapped_rejected += eval_stats.gapped_rejected;

        for decision in passes {
            let sid = decision.seq_id;
            if decision.passed && assignments[sid as usize].is_none() {
                assignments[sid as usize] = Some(cid);
                counters.assigned_pairs += 1;
            }
        }
    }
    drop(greedy_guard);

    let assignments = assignments
        .into_iter()
        .enumerate()
        .map(|(i, cid)| {
            cid.unwrap_or_else(|| {
                let new_cid = representatives.len() as u32;
                representatives.push(i as u32);
                new_cid
            })
        })
        .collect();

    if let Some(p) = &cfg.profiler {
        p.add_counters(&counters);
    }

    Ok(ClusterResult {
        representatives,
        assignments,
    })
}

fn build_group_candidates(db: &SeqDb, table: &KmerTable) -> Vec<Vec<u32>> {
    let mut out: Vec<Vec<u32>> = vec![Vec::new(); db.n_seqs()];

    for group in table.groups() {
        if group.len() < 2 {
            continue;
        }

        let mut ids: Vec<u32> = group.iter().map(|e| e.seq_id).collect();
        ids.sort_unstable();
        ids.dedup();
        if ids.len() < 2 {
            continue;
        }

        ids.sort_by(|a, b| db.seq_len(*b).cmp(&db.seq_len(*a)).then(a.cmp(b)));
        let leader = ids[0];
        out[leader as usize].extend_from_slice(&ids[1..]);
    }

    for v in &mut out {
        v.sort_unstable();
        v.dedup();
    }

    out
}
