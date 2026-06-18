use std::fs;
use std::path::Path;

use crate::error::{AppError, Result};
use crate::index::build::{BuildIndexConfig, build_next_segment};
use crate::index::manifest::{ManifestSegment, load_manifest, save_manifest};
use crate::index::segment::{SegmentMeta, segment_dir};

pub fn update_index(indexdir: &Path, cfg: BuildIndexConfig) -> Result<u32> {
    fs::create_dir_all(indexdir)
        .map_err(|e| AppError::io(format!("create dir {}", indexdir.display()), e))?;

    let mut manifest = load_manifest(indexdir)?;
    let global_seq_base = next_global_seq_base(indexdir, &manifest.segments)?;
    let seg_id = build_next_segment(&cfg, &manifest, global_seq_base)?;

    manifest
        .segments
        .push(ManifestSegment { segment_id: seg_id });
    manifest.segments.sort_by_key(|s| s.segment_id);
    save_manifest(indexdir, &manifest)?;

    Ok(seg_id)
}

fn next_global_seq_base(indexdir: &Path, segments: &[ManifestSegment]) -> Result<u64> {
    let mut next = 0u64;
    for seg in segments {
        let path = segment_dir(indexdir, seg.segment_id).join("segment_meta.json");
        let bytes = std::fs::read(&path)
            .map_err(|e| AppError::io(format!("read {}", path.display()), e))?;
        let meta: SegmentMeta = serde_json::from_slice(&bytes).map_err(|e| {
            AppError::new(
                crate::error::ErrorKind::Parse,
                format!("parse {}: {e}", path.display()),
            )
        })?;
        next = next.max(meta.global_seq_base.saturating_add(meta.n_seqs));
    }
    Ok(next)
}
