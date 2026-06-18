use std::path::{Path, PathBuf};

pub fn is_db_dir(path: &Path) -> bool {
    path.is_dir() && path.join("meta.json").exists() && path.join("seqs.bin").exists()
}

pub fn query_db_dir(tmp_dir: &Path) -> PathBuf {
    tmp_dir.join("query_db")
}
