#[derive(Debug, Clone)]
pub struct ScoringParams {
    pub match_score: i16,
    pub mismatch_score: i16,
    pub gap_open: i16,
    pub gap_extend: i16,
    /// When true, `score_pair` uses the BLOSUM62 substitution matrix instead of
    /// the flat `match_score`/`mismatch_score` model. This is the default for
    /// protein scoring and matches the behaviour of MMseqs2/BLAST/CD-HIT.
    pub use_matrix: bool,
}

impl ScoringParams {
    /// Flat match/mismatch protein scoring. Kept for nucleotide-like or
    /// substitution-agnostic use and for low-level alignment unit tests.
    pub fn protein_default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -1,
            gap_open: -11,
            gap_extend: -1,
            use_matrix: false,
        }
    }

    /// Protein scoring backed by BLOSUM62 with standard affine gap penalties.
    pub fn protein_blosum62() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -1,
            gap_open: -11,
            gap_extend: -1,
            use_matrix: true,
        }
    }

    pub fn nucleotide_default() -> Self {
        Self {
            match_score: 2,
            mismatch_score: -2,
            gap_open: -8,
            gap_extend: -2,
            use_matrix: false,
        }
    }

    pub fn score_pair(&self, a: u8, b: u8) -> i16 {
        if self.use_matrix {
            blosum62(a, b)
        } else if a == b {
            self.match_score
        } else {
            self.mismatch_score
        }
    }
}

/// BLOSUM62 symbol order (NCBI standard, excluding the `*` stop column).
const BLOSUM62_ORDER: &[u8; 23] = b"ARNDCQEGHILKMFPSTWYVBZX";
const X_INDEX: u8 = 22;

#[rustfmt::skip]
const BLOSUM62_MATRIX: [[i8; 23]; 23] = [
    //   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X
    [    4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0], // A
    [   -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1], // R
    [   -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1], // N
    [   -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1], // D
    [    0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2], // C
    [   -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1], // Q
    [   -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1], // E
    [    0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1], // G
    [   -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1], // H
    [   -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1], // I
    [   -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1], // L
    [   -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1], // K
    [   -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1], // M
    [   -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1], // F
    [   -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2], // P
    [    1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0], // S
    [    0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0], // T
    [   -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2], // W
    [   -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1], // Y
    [    0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1], // V
    [   -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1], // B
    [   -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1], // Z
    [    0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1], // X
];

/// Maps an ASCII residue byte to its BLOSUM62 row/column index. Unknown bytes
/// (including `*`, gaps, digits) fall back to the `X` column.
const fn build_blosum_index() -> [u8; 256] {
    let mut idx = [X_INDEX; 256];
    let mut i = 0;
    while i < BLOSUM62_ORDER.len() {
        let c = BLOSUM62_ORDER[i];
        idx[c as usize] = i as u8;
        idx[(c + 32) as usize] = i as u8; // lowercase
        i += 1;
    }
    idx
}

const BLOSUM_INDEX: [u8; 256] = build_blosum_index();

#[inline]
pub fn blosum62(a: u8, b: u8) -> i16 {
    let ia = BLOSUM_INDEX[a as usize] as usize;
    let ib = BLOSUM_INDEX[b as usize] as usize;
    BLOSUM62_MATRIX[ia][ib] as i16
}

#[cfg(test)]
mod tests {
    use super::{ScoringParams, blosum62};

    #[test]
    fn blosum62_is_symmetric_and_known_values() {
        assert_eq!(blosum62(b'A', b'A'), 4);
        assert_eq!(blosum62(b'W', b'W'), 11);
        assert_eq!(blosum62(b'C', b'C'), 9);
        // Conservative substitution scores positive.
        assert_eq!(blosum62(b'Q', b'R'), 1);
        assert_eq!(blosum62(b'I', b'L'), 2);
        // Distant substitution scores negative.
        assert_eq!(blosum62(b'W', b'P'), -4);
        // Symmetry across the whole alphabet.
        for &x in super::BLOSUM62_ORDER {
            for &y in super::BLOSUM62_ORDER {
                assert_eq!(blosum62(x, y), blosum62(y, x));
            }
        }
    }

    #[test]
    fn unknown_bytes_fall_back_to_x() {
        assert_eq!(blosum62(b'*', b'A'), blosum62(b'X', b'A'));
        assert_eq!(blosum62(b'A', b'a'), blosum62(b'A', b'A'));
    }

    #[test]
    fn score_pair_respects_use_matrix() {
        let flat = ScoringParams::protein_default();
        assert_eq!(flat.score_pair(b'Q', b'R'), -1);
        let blosum = ScoringParams::protein_blosum62();
        assert_eq!(blosum.score_pair(b'Q', b'R'), 1);
    }
}
