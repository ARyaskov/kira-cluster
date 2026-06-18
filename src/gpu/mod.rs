use crate::alignment::scoring::ScoringParams;
use crate::alignment::{AlignmentResult, SeqId};
use crate::cluster::SeqDb;
use crate::error::{AppError, ErrorKind, Result};

pub mod cuda;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GpuBackend {
    Cuda,
}

impl GpuBackend {
    pub fn as_str(self) -> &'static str {
        match self {
            GpuBackend::Cuda => "cuda",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct GpuConfig {
    pub backend: GpuBackend,
    pub batch_size: usize,
    pub gpu_memory_limit: usize,
}

pub fn is_backend_available(backend: GpuBackend) -> bool {
    match backend {
        GpuBackend::Cuda => cuda::is_available(),
    }
}

pub fn align_pairs_on_gpu(
    db: &SeqDb,
    pairs: &[(SeqId, SeqId)],
    scoring: &ScoringParams,
    cfg: GpuConfig,
) -> Result<Vec<AlignmentResult>> {
    let _ = cfg.batch_size;
    let _ = cfg.gpu_memory_limit;

    if pairs.is_empty() {
        return Ok(Vec::new());
    }

    match cfg.backend {
        GpuBackend::Cuda => {
            if !cuda::is_available() {
                return Err(AppError::new(
                    ErrorKind::Unsupported,
                    "CUDA backend unavailable",
                ));
            }
            cuda::align_pairs_cuda(db, pairs, scoring)
        }
    }
}
