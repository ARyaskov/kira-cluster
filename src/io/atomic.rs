use std::fs;
use std::path::{Path, PathBuf};

use crate::error::{AppError, Result};

fn tmp_path_for(path: &Path) -> PathBuf {
    let fname = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("tmpfile");
    path.with_file_name(format!("{fname}.tmp"))
}

pub fn write_atomic(path: &Path, bytes: &[u8]) -> Result<()> {
    let tmp = tmp_path_for(path);
    fs::write(&tmp, bytes).map_err(|e| AppError::io(format!("write {}", tmp.display()), e))?;
    fs::rename(&tmp, path)
        .map_err(|e| AppError::io(format!("rename {} -> {}", tmp.display(), path.display()), e))?;
    Ok(())
}
