use crate::alignment::AlignmentResult;
use crate::alignment::gapped_simd::smith_waterman_gotoh_simd_neon;
use crate::alignment::scoring::ScoringParams;

pub fn smith_waterman_gotoh_neon(a: &[u8], b: &[u8], scoring: &ScoringParams) -> AlignmentResult {
    smith_waterman_gotoh_simd_neon(a, b, scoring)
}
