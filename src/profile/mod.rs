use std::collections::BTreeMap;
use std::path::Path;
use std::sync::Mutex;
use std::time::Instant;

use serde::Serialize;

use crate::error::{AppError, ErrorKind, Result};
use crate::io::atomic::write_atomic;

pub mod timer;

#[derive(Debug, Clone, Default, Serialize)]
pub struct WorkCounters {
    pub n_sequences: u64,
    pub total_seeds_emitted: u64,
    pub unique_keys: u64,
    pub total_posting_length_read: u64,
    pub candidate_pairs_evaluated: u64,
    pub hamming_rejected: u64,
    pub ungapped_rejected: u64,
    pub gapped_rejected: u64,
    pub assigned_pairs: u64,
}

impl WorkCounters {
    pub fn merge(&mut self, rhs: &WorkCounters) {
        self.n_sequences += rhs.n_sequences;
        self.total_seeds_emitted += rhs.total_seeds_emitted;
        self.unique_keys += rhs.unique_keys;
        self.total_posting_length_read += rhs.total_posting_length_read;
        self.candidate_pairs_evaluated += rhs.candidate_pairs_evaluated;
        self.hamming_rejected += rhs.hamming_rejected;
        self.ungapped_rejected += rhs.ungapped_rejected;
        self.gapped_rejected += rhs.gapped_rejected;
        self.assigned_pairs += rhs.assigned_pairs;
    }
}

#[derive(Debug)]
struct ProfileState {
    stage_timings_ns: BTreeMap<&'static str, u128>,
    counters: WorkCounters,
}

impl Default for ProfileState {
    fn default() -> Self {
        Self {
            stage_timings_ns: BTreeMap::new(),
            counters: WorkCounters::default(),
        }
    }
}

#[derive(Debug)]
pub struct Profiler {
    enabled: bool,
    state: Mutex<ProfileState>,
    wall_start: Instant,
}

impl Profiler {
    pub fn new(enabled: bool) -> Self {
        Self {
            enabled,
            state: Mutex::new(ProfileState::default()),
            wall_start: Instant::now(),
        }
    }

    pub fn enabled(&self) -> bool {
        self.enabled
    }

    pub fn stage(&self, name: &'static str) -> StageGuard<'_> {
        if !self.enabled {
            return StageGuard {
                profiler: None,
                name,
                start: Instant::now(),
            };
        }
        StageGuard {
            profiler: Some(self),
            name,
            start: Instant::now(),
        }
    }

    pub fn add_stage_time(&self, name: &'static str, elapsed_ns: u128) {
        if !self.enabled {
            return;
        }
        let mut state = self.state.lock().expect("profile mutex poisoned");
        *state.stage_timings_ns.entry(name).or_insert(0) += elapsed_ns;
    }

    pub fn add_counters(&self, counters: &WorkCounters) {
        if !self.enabled {
            return;
        }
        let mut state = self.state.lock().expect("profile mutex poisoned");
        state.counters.merge(counters);
    }

    pub fn set_unique_keys(&self, n: u64) {
        if !self.enabled {
            return;
        }
        let mut state = self.state.lock().expect("profile mutex poisoned");
        state.counters.unique_keys = n;
    }

    pub fn set_n_sequences(&self, n: u64) {
        if !self.enabled {
            return;
        }
        let mut state = self.state.lock().expect("profile mutex poisoned");
        state.counters.n_sequences = n;
    }

    pub fn snapshot(&self) -> ProfileReport {
        let mut state = self.state.lock().expect("profile mutex poisoned");
        let total_wall_clock_ns = self.wall_start.elapsed().as_nanos();
        state.stage_timings_ns.entry("parse_ingest").or_insert(0);
        state.stage_timings_ns.entry("seed_generation").or_insert(0);
        state
            .stage_timings_ns
            .entry("kmer_table_build")
            .or_insert(0);
        state.stage_timings_ns.entry("kmer_sort").or_insert(0);
        state
            .stage_timings_ns
            .entry("candidate_grouping")
            .or_insert(0);
        state
            .stage_timings_ns
            .entry("filter_prefilter")
            .or_insert(0);
        state.stage_timings_ns.entry("filter_gapped").or_insert(0);
        state
            .stage_timings_ns
            .entry("greedy_assignment")
            .or_insert(0);
        state.stage_timings_ns.entry("final_write").or_insert(0);

        let c = state.counters.clone();
        let surv_after_hamming = c
            .candidate_pairs_evaluated
            .saturating_sub(c.hamming_rejected);
        let surv_after_ungapped = surv_after_hamming.saturating_sub(c.ungapped_rejected);

        ProfileReport {
            n_sequences: c.n_sequences,
            total_seeds_emitted: c.total_seeds_emitted,
            unique_keys: c.unique_keys,
            total_posting_length_read: c.total_posting_length_read,
            candidate_pairs_evaluated: c.candidate_pairs_evaluated,
            hamming_rejected: c.hamming_rejected,
            ungapped_rejected: c.ungapped_rejected,
            gapped_rejected: c.gapped_rejected,
            assigned_pairs: c.assigned_pairs,
            reject_rates: RejectRates {
                hamming: ratio(c.hamming_rejected, c.candidate_pairs_evaluated),
                ungapped: ratio(c.ungapped_rejected, surv_after_hamming),
                gapped: ratio(c.gapped_rejected, surv_after_ungapped),
            },
            stage_timings_ns: state
                .stage_timings_ns
                .iter()
                .map(|(k, v)| ((*k).to_string(), *v))
                .collect(),
            total_wall_clock_ns,
        }
    }
}

pub struct StageGuard<'a> {
    profiler: Option<&'a Profiler>,
    name: &'static str,
    start: Instant,
}

impl Drop for StageGuard<'_> {
    fn drop(&mut self) {
        if let Some(p) = self.profiler {
            p.add_stage_time(self.name, self.start.elapsed().as_nanos());
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct RejectRates {
    pub hamming: f64,
    pub ungapped: f64,
    pub gapped: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ProfileReport {
    pub n_sequences: u64,
    pub total_seeds_emitted: u64,
    pub unique_keys: u64,
    pub total_posting_length_read: u64,
    pub candidate_pairs_evaluated: u64,
    pub hamming_rejected: u64,
    pub ungapped_rejected: u64,
    pub gapped_rejected: u64,
    pub assigned_pairs: u64,
    pub reject_rates: RejectRates,
    pub stage_timings_ns: BTreeMap<String, u128>,
    pub total_wall_clock_ns: u128,
}

impl ProfileReport {
    pub fn print_console(&self) {
        println!("=== kira-cluster profiling ===");
        println!("N sequences: {}", self.n_sequences);
        println!("Total seeds emitted: {}", self.total_seeds_emitted);
        println!("Unique keys: {}", self.unique_keys);
        println!(
            "Candidate pairs evaluated: {}",
            self.candidate_pairs_evaluated
        );
        println!();
        println!("Reject rates:");
        println!("  Hamming: {:.3}%", self.reject_rates.hamming * 100.0);
        println!("  Ungapped: {:.3}%", self.reject_rates.ungapped * 100.0);
        println!("  Gapped: {:.3}%", self.reject_rates.gapped * 100.0);
        println!();
        println!("Stage timings (ms):");
        for name in [
            "parse_ingest",
            "seed_generation",
            "kmer_table_build",
            "kmer_sort",
            "candidate_grouping",
            "filter_prefilter",
            "filter_gapped",
            "greedy_assignment",
            "final_write",
        ] {
            let ms = self.stage_timings_ns.get(name).copied().unwrap_or(0) / 1_000_000;
            println!("  {:<17} {}", name, ms);
        }
        println!(
            "Total wall-clock: {} ms",
            self.total_wall_clock_ns / 1_000_000
        );
    }

    pub fn write_json(&self, path: &Path) -> Result<()> {
        let bytes = serde_json::to_vec_pretty(self).map_err(|e| {
            AppError::new(ErrorKind::Internal, format!("serialize profile json: {e}"))
        })?;
        write_atomic(path, &bytes)
    }
}

fn ratio(num: u64, den: u64) -> f64 {
    if den == 0 {
        0.0
    } else {
        (num as f64) / (den as f64)
    }
}
