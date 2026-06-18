use crate::alignment::AlignmentResult;
use crate::alignment::gapped_simd::smith_waterman_gotoh_simd_avx2;
use crate::alignment::scoring::ScoringParams;

pub fn smith_waterman_gotoh_avx2(a: &[u8], b: &[u8], scoring: &ScoringParams) -> AlignmentResult {
    smith_waterman_gotoh_simd_avx2(a, b, scoring)
}
