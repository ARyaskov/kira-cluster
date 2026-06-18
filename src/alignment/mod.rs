use crate::alignment::scoring::ScoringParams;
use crate::simd::SimdBackend;

pub mod gapped_avx2;
pub mod gapped_neon;
pub mod gapped_scalar;
mod gapped_simd;
pub mod gpu;
pub mod scoring;

pub type SeqId = u32;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentMode {
    Fast,
    Sensitive,
}

impl AlignmentMode {
    pub fn as_str(self) -> &'static str {
        match self {
            AlignmentMode::Fast => "fast",
            AlignmentMode::Sensitive => "sensitive",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignmentResult {
    pub score: i32,
    pub aligned_len: u32,
    pub matches: u32,
}

pub fn smith_waterman_gotoh(
    a: &[u8],
    b: &[u8],
    scoring: &ScoringParams,
    backend: SimdBackend,
) -> AlignmentResult {
    match backend {
        SimdBackend::Avx2 => gapped_avx2::smith_waterman_gotoh_avx2(a, b, scoring),
        SimdBackend::Neon => gapped_neon::smith_waterman_gotoh_neon(a, b, scoring),
        SimdBackend::Scalar => gapped_scalar::smith_waterman_gotoh_scalar(a, b, scoring),
    }
}
