use crate::alignment::{AlignmentResult, SeqId, scoring::ScoringParams};
use crate::error::{AppError, ErrorKind, Result};

pub trait DbView {
    fn seq(&self, seq_id: SeqId) -> &[u8];
}

pub fn smith_waterman_gotoh_gpu(
    _pairs: &[(SeqId, SeqId)],
    _db: &dyn DbView,
    _scoring: &ScoringParams,
) -> Result<Vec<AlignmentResult>> {
    Err(AppError::new(
        ErrorKind::Unsupported,
        "GPU backend not yet available",
    ))
}
