use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::error::{AppError, Result};

fn tmp_path_for(path: &Path) -> PathBuf {
    let fname = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("tmpfile");
    path.with_file_name(format!("{fname}.tmp"))
}

/// Ensure the parent directory of `path` exists.
///
/// Atomic writes create a sibling `*.tmp` file next to the destination, so the
/// destination's parent directory must exist first. Callers such as
/// `--profile-json` previously failed with `KC1002 ... No such file or directory`
/// when the parent directory had not been created.
fn ensure_parent_dir(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .map_err(|e| AppError::io(format!("create dir {}", parent.display()), e))?;
        }
    }
    Ok(())
}

pub fn write_atomic(path: &Path, bytes: &[u8]) -> Result<()> {
    ensure_parent_dir(path)?;
    let tmp = tmp_path_for(path);
    fs::write(&tmp, bytes).map_err(|e| AppError::io(format!("write {}", tmp.display()), e))?;
    fs::rename(&tmp, path)
        .map_err(|e| AppError::io(format!("rename {} -> {}", tmp.display(), path.display()), e))?;
    Ok(())
}

/// Stream bytes into `path` atomically without materializing the whole payload
/// in memory. The closure receives a buffered writer over a temporary file that
/// is renamed into place on success.
pub fn write_atomic_with<F>(path: &Path, write_body: F) -> Result<()>
where
    F: FnOnce(&mut BufWriter<File>) -> Result<()>,
{
    ensure_parent_dir(path)?;
    let tmp = tmp_path_for(path);
    {
        let file =
            File::create(&tmp).map_err(|e| AppError::io(format!("create {}", tmp.display()), e))?;
        let mut writer = BufWriter::new(file);
        write_body(&mut writer)?;
        writer
            .flush()
            .map_err(|e| AppError::io(format!("flush {}", tmp.display()), e))?;
    }
    fs::rename(&tmp, path)
        .map_err(|e| AppError::io(format!("rename {} -> {}", tmp.display(), path.display()), e))?;
    Ok(())
}
