#[derive(Debug, Clone)]
pub struct ScoringParams {
    pub match_score: i16,
    pub mismatch_score: i16,
    pub gap_open: i16,
    pub gap_extend: i16,
}

impl ScoringParams {
    pub fn protein_default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -1,
            gap_open: -11,
            gap_extend: -1,
        }
    }

    pub fn nucleotide_default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -2,
            gap_open: -8,
            gap_extend: -2,
        }
    }

    pub fn score_pair(&self, a: u8, b: u8) -> i16 {
        if a == b {
            self.match_score
        } else {
            self.mismatch_score
        }
    }
}

pub fn score(a: u8, b: u8) -> i16 {
    if a == b { 2 } else { -1 }
}
