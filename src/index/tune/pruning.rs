use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::{AppError, ErrorKind, Result};
use crate::io::atomic::write_atomic;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SegmentStats {
    pub version: u32,
    pub segment_id: u32,
    pub n_seqs: u64,
    pub n_keys: u64,
    pub min_key: u64,
    pub max_key: u64,
    pub df_p50: u32,
    pub df_p90: u32,
    pub df_max: u32,
}

impl SegmentStats {
    pub fn from_data(segment_id: u32, n_seqs: u64, keys: &[u64], dfs: &[u32]) -> Self {
        let mut sorted = dfs.to_vec();
        sorted.sort_unstable();
        let q = |p: f32| -> u32 {
            if sorted.is_empty() {
                return 0;
            }
            let idx = ((sorted.len() - 1) as f32 * p).round() as usize;
            sorted[idx.min(sorted.len() - 1)]
        };

        Self {
            version: 1,
            segment_id,
            n_seqs,
            n_keys: keys.len() as u64,
            min_key: keys.first().copied().unwrap_or(0),
            max_key: keys.last().copied().unwrap_or(0),
            df_p50: q(0.50),
            df_p90: q(0.90),
            df_max: *sorted.last().unwrap_or(&0),
        }
    }
}

pub fn write_segment_stats(segment_dir: &Path, stats: &SegmentStats) -> Result<()> {
    let bytes = serde_json::to_vec_pretty(stats)
        .map_err(|e| AppError::new(ErrorKind::Internal, format!("serialize segment stats: {e}")))?;
    write_atomic(&segment_dir.join("segment_stats.json"), &bytes)
}
