use crate::alignment::SeqId;
use crate::alignment::scoring::ScoringParams;

#[derive(Debug, Clone)]
pub struct GpuBatch {
    pub pairs: Vec<(SeqId, SeqId)>,
    pub scoring: ScoringParams,
}

pub fn build_batches(
    pairs: &[(SeqId, SeqId)],
    scoring: &ScoringParams,
    batch_size: usize,
) -> Vec<GpuBatch> {
    let size = batch_size.max(1);
    let mut out = Vec::new();
    for chunk in pairs.chunks(size) {
        out.push(GpuBatch {
            pairs: chunk.to_vec(),
            scoring: scoring.clone(),
        });
    }
    out
}
