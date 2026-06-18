use crate::alignment::scoring::ScoringParams;
use crate::alignment::smith_waterman_gotoh;
use crate::cluster::filter::{hamming_filter, ungapped_filter};
use crate::simd::SimdBackend;

#[derive(Debug, Clone, Copy)]
pub struct CandidateScore {
    pub score: i32,
    pub identity: f32,
    pub coverage: f32,
}

pub fn score_candidate(
    query: &[u8],
    target: &[u8],
    seeds_hit: u32,
    min_identity: f32,
    min_coverage: f32,
    scoring: &ScoringParams,
    backend: SimdBackend,
) -> Option<CandidateScore> {
    let overlap = query.len().min(target.len()) as f32;
    let allowed = ((1.0 - min_identity) * overlap).floor().max(0.0) as u32;

    if !hamming_filter(query, target, allowed + 16, backend) {
        return None;
    }
    if !ungapped_filter(
        query,
        target,
        (min_identity * 0.9).clamp(0.0, 1.0),
        (min_coverage * 0.8).clamp(0.0, 1.0),
        backend,
    ) {
        return None;
    }

    let aln = smith_waterman_gotoh(query, target, scoring, backend);
    if aln.aligned_len == 0 {
        return None;
    }

    let identity = (aln.matches as f32) / (aln.aligned_len as f32);
    let coverage = (aln.aligned_len as f32) / (query.len().min(target.len()) as f32);
    if identity < min_identity || coverage < min_coverage {
        return None;
    }

    Some(CandidateScore {
        score: aln.score + (seeds_hit as i32) * 4,
        identity,
        coverage,
    })
}
