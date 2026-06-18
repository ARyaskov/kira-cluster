use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use crate::error::{AppError, ErrorKind, Result};
use crate::index::learned::PgmIndexV2;
use crate::index::manifest::load_manifest;
use crate::index::segment::{SegmentHandle, segment_dir};
use crate::index::tune::pruning::{SegmentStats, write_segment_stats};
use crate::index::tune::seed_policy::{SeedPolicy, choose_default_policy, write_seed_policy};

pub mod pruning;
pub mod seed_policy;

#[derive(Debug, Clone)]
pub struct TuneConfig {
    pub indexdir: PathBuf,
    pub pgm_epsilon: usize,
    pub seed_policy_mode: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TuneMeta {
    pub version: u32,
    pub pgm_epsilon: u32,
    pub tuned_segments: u32,
}

pub fn tune_index(cfg: &TuneConfig) -> Result<()> {
    let manifest = load_manifest(&cfg.indexdir)?;

    for seg in &manifest.segments {
        let dir = segment_dir(&cfg.indexdir, seg.segment_id);
        let handle = SegmentHandle::open(&cfg.indexdir, seg.segment_id)?;
        let keys = handle.keys();

        let pgm = PgmIndexV2::build(keys, cfg.pgm_epsilon)?;
        pgm.save_to_segment(&dir)?;

        let dfs = collect_dfs(&handle, keys.len());
        let stats = SegmentStats::from_data(handle.meta.segment_id, handle.meta.n_seqs, keys, &dfs);
        write_segment_stats(&dir, &stats)?;
    }

    let policy = if let Some(mode) = &cfg.seed_policy_mode {
        SeedPolicy::from_mode_str(mode).ok_or_else(|| {
            AppError::new(
                ErrorKind::Validation,
                format!("invalid seed policy mode: {mode}"),
            )
        })?
    } else {
        choose_default_policy()
    };
    write_seed_policy(&cfg.indexdir, &policy)?;

    let meta = TuneMeta {
        version: 1,
        pgm_epsilon: cfg.pgm_epsilon.max(1) as u32,
        tuned_segments: manifest.segments.len() as u32,
    };
    let bytes = serde_json::to_vec_pretty(&meta)
        .map_err(|e| AppError::new(ErrorKind::Internal, format!("serialize tune meta: {e}")))?;
    crate::io::atomic::write_atomic(&cfg.indexdir.join("tune_meta.json"), &bytes)
}

fn collect_dfs(handle: &SegmentHandle, n_keys: usize) -> Vec<u32> {
    let mut out = Vec::with_capacity(n_keys);
    for i in 0..n_keys {
        out.push(handle.df(i));
    }
    out
}
