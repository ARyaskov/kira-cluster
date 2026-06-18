use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::error::{AppError, ErrorKind, Result};
use crate::io::atomic::write_atomic;
use crate::io::mmap::map_readonly;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ManifestSegment {
    pub segment_id: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Manifest {
    pub version: u32,
    pub index_name: String,
    pub segments: Vec<ManifestSegment>,
}

pub fn manifest_path(indexdir: &Path) -> PathBuf {
    indexdir.join("manifest.json")
}

pub fn load_manifest(indexdir: &Path) -> Result<Manifest> {
    let path = manifest_path(indexdir);
    let map = map_readonly(&path)?;
    serde_json::from_slice(&map)
        .map_err(|e| AppError::new(ErrorKind::Parse, format!("parse {}: {e}", path.display())))
}

pub fn save_manifest(indexdir: &Path, manifest: &Manifest) -> Result<()> {
    let bytes = serde_json::to_vec_pretty(manifest)
        .map_err(|e| AppError::new(ErrorKind::Internal, format!("serialize manifest: {e}")))?;
    write_atomic(&manifest_path(indexdir), &bytes)
}

pub fn next_segment_id(manifest: &Manifest) -> u32 {
    manifest
        .segments
        .iter()
        .map(|s| s.segment_id)
        .max()
        .map_or(0, |m| m + 1)
}
