use std::path::Path;

use kira_kv_engine::PgmIndex;

use crate::error::{AppError, ErrorKind, Result};
use crate::io::atomic::write_atomic;
use crate::io::mmap::map_readonly;

const KIRA_KV_PGM_FILE: &str = "kv_pgm.bin";

pub struct PgmIndexV2 {
    inner: PgmIndex,
}

impl PgmIndexV2 {
    pub fn build(sorted_keys: &[u64], epsilon: usize) -> Result<Self> {
        let epsilon = normalize_epsilon(epsilon);
        let inner = PgmIndex::build(sorted_keys.to_vec(), epsilon).map_err(|e| {
            AppError::new(
                ErrorKind::Internal,
                format!("build kira_kv_engine PGM: {e}"),
            )
        })?;
        Ok(Self { inner })
    }

    pub fn find_key(&self, key: u64) -> Option<usize> {
        self.inner.index(key).ok()
    }

    pub fn save_to_segment(&self, segment_dir: &Path) -> Result<()> {
        let bytes = self.inner.to_bytes().map_err(|e| {
            AppError::new(
                ErrorKind::Internal,
                format!(
                    "serialize kira_kv_engine PGM in {}: {e}",
                    segment_dir.display()
                ),
            )
        })?;
        write_atomic(&segment_dir.join(KIRA_KV_PGM_FILE), &bytes)
    }

    pub fn load_from_segment(segment_dir: &Path) -> Result<Self> {
        let bytes = map_readonly(&segment_dir.join(KIRA_KV_PGM_FILE))?;
        let inner = PgmIndex::from_bytes(&bytes).map_err(|e| {
            AppError::new(ErrorKind::Parse, format!("parse kira_kv_engine PGM: {e}"))
        })?;
        Ok(Self { inner })
    }
}

fn normalize_epsilon(epsilon: usize) -> u32 {
    u32::try_from(epsilon.max(1)).unwrap_or(u32::MAX)
}
