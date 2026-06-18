use std::fs::File;
use std::path::Path;

use memmap2::Mmap;

use crate::error::{AppError, Result};

pub fn map_readonly(path: &Path) -> Result<Mmap> {
    let file = File::open(path).map_err(|e| AppError::io(format!("open {}", path.display()), e))?;
    // SAFETY: read-only mapping of an immutable file descriptor in this process.
    let map = unsafe { Mmap::map(&file) }
        .map_err(|e| AppError::io(format!("mmap {}", path.display()), e))?;
    Ok(map)
}
