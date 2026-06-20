//! Banded local Gotoh alignment with a score-only forward pass and traceback.
//!
//! The forward pass computes only the H/E/F **scores** (plain i32 max/add, no
//! per-cell `aligned_len`/`matches` bookkeeping) and records a one-byte
//! traceback direction per banded cell. `matches`/`aligned_len` — needed for the
//! identity/coverage gate — are then recovered by a single O(alignment length)
//! traceback. This is much cheaper per cell than carrying a `Cell { score,
//! aligned_len, matches }` through the recurrence, and the score-only forward is
//! the part a striped-SIMD (Farrar) kernel would later replace.

use crate::alignment::AlignmentResult;
use crate::alignment::scoring::ScoringParams;

const NEG_INF: i32 = i32::MIN / 4;

// Traceback direction encoding (one byte per banded cell).
const H_ZERO: u8 = 0; // local restart (cell is the start boundary)
const H_DIAG: u8 = 1; // H came from the diagonal (a match/mismatch column)
const H_FROM_E: u8 = 2; // H came from the E (gap-in-a) state
const H_FROM_F: u8 = 3; // H came from the F (gap-in-b) state
const H_MASK: u8 = 0b11;
const E_EXTEND: u8 = 1 << 2; // else E was opened from H
const F_EXTEND: u8 = 1 << 3; // else F was opened from H

/// Largest direction matrix (rows * band columns) we will allocate. Beyond this
/// the pair is treated as un-alignable (returned as a zero alignment). The
/// coverage prefilter keeps near-equal lengths, so this is only a safety valve.
const MAX_DIR_CELLS: usize = 64_000_000;

pub fn smith_waterman_gotoh_banded_tb(
    a: &[u8],
    b: &[u8],
    scoring: &ScoringParams,
    half_band: usize,
) -> AlignmentResult {
    if a.is_empty() || b.is_empty() {
        return zero();
    }

    let m = a.len();
    let n = b.len();
    let go_ge = scoring.gap_open as i32 + scoring.gap_extend as i32;
    let ge = scoring.gap_extend as i32;

    // Static band on the column offset d = j - i (see gapped_scalar).
    let w = half_band as isize;
    let d_end = n as isize - m as isize;
    let d_lo = d_end.min(0) - w;
    let d_hi = d_end.max(0) + w;
    let band_span = (d_hi - d_lo + 1) as usize;

    let dir_cells = match m.checked_mul(band_span) {
        Some(c) if c <= MAX_DIR_CELLS => c,
        _ => return zero(),
    };
    let mut dir = vec![0u8; dir_cells];

    let mut prev_h = vec![0i32; n + 1];
    let mut curr_h = vec![0i32; n + 1];
    let mut prev_f = vec![NEG_INF; n + 1];
    let mut curr_f = vec![NEG_INF; n + 1];

    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;
    let mut prev_hi = 0usize;

    for i in 1..=m {
        let ii = i as isize;
        let lo = (ii + d_lo).max(1) as usize;
        let hi_signed = (ii + d_hi).min(n as isize);
        if hi_signed < lo as isize {
            for x in curr_h.iter_mut() {
                *x = 0;
            }
            for x in curr_f.iter_mut() {
                *x = NEG_INF;
            }
            std::mem::swap(&mut prev_h, &mut curr_h);
            std::mem::swap(&mut prev_f, &mut curr_f);
            prev_hi = 0;
            continue;
        }
        let hi = hi_signed as usize;

        for jj in (prev_hi + 1)..=hi {
            prev_h[jj] = 0;
            prev_f[jj] = NEG_INF;
        }
        curr_h[lo - 1] = 0;

        let row_base = (i - 1) * band_span;
        let col0 = ii + d_lo; // band index of column j is (j - col0)
        let ai = a[i - 1];
        let mut e = NEG_INF;
        for j in lo..=hi {
            let sub = scoring.score_pair(ai, b[j - 1]) as i32;
            let diag = prev_h[j - 1] + sub; // prev_h >= 0 (local)

            // E: horizontal gap (consumes b).
            let e_open = curr_h[j - 1] + go_ge;
            let e_ext = e + ge;
            let e_bit;
            if e_ext > e_open {
                e = e_ext;
                e_bit = E_EXTEND;
            } else {
                e = e_open;
                e_bit = 0;
            }

            // F: vertical gap (consumes a).
            let f_open = prev_h[j] + go_ge;
            let f_ext = prev_f[j] + ge;
            let f_bit;
            let f_val;
            if f_ext > f_open {
                f_val = f_ext;
                f_bit = F_EXTEND;
            } else {
                f_val = f_open;
                f_bit = 0;
            }
            curr_f[j] = f_val;

            // H = max(0, diag, E, F). Deterministic preference diag > E > F.
            let mut h = 0i32;
            let mut hsrc = H_ZERO;
            if diag > h {
                h = diag;
                hsrc = H_DIAG;
            }
            if e > h {
                h = e;
                hsrc = H_FROM_E;
            }
            if f_val > h {
                h = f_val;
                hsrc = H_FROM_F;
            }
            curr_h[j] = h;

            let band_idx = (j as isize - col0) as usize;
            dir[row_base + band_idx] = hsrc | e_bit | f_bit;

            if h > best_score {
                best_score = h;
                best_i = i;
                best_j = j;
            }
        }

        prev_hi = hi;
        std::mem::swap(&mut prev_h, &mut curr_h);
        std::mem::swap(&mut prev_f, &mut curr_f);
    }

    if best_score <= 0 {
        return zero();
    }

    // Traceback for matches / aligned_len.
    let mut matches = 0u32;
    let mut aligned_len = 0u32;
    let mut i = best_i;
    let mut j = best_j;
    let mut state = State::H;

    loop {
        if i == 0 || j == 0 {
            break;
        }
        let col0 = i as isize + d_lo;
        let band_idx = j as isize - col0;
        if band_idx < 0 || band_idx as usize >= band_span {
            break; // out of band: should not happen for an in-band optimum
        }
        let d = dir[(i - 1) * band_span + band_idx as usize];

        match state {
            State::H => match d & H_MASK {
                H_ZERO => break,
                H_DIAG => {
                    aligned_len += 1;
                    if a[i - 1] == b[j - 1] {
                        matches += 1;
                    }
                    i -= 1;
                    j -= 1;
                }
                H_FROM_E => state = State::E,
                _ => state = State::F,
            },
            State::E => {
                aligned_len += 1; // gap in a, consumes a b column
                let ext = d & E_EXTEND != 0;
                j -= 1;
                if !ext {
                    state = State::H;
                }
            }
            State::F => {
                aligned_len += 1; // gap in b, consumes an a column
                let ext = d & F_EXTEND != 0;
                i -= 1;
                if !ext {
                    state = State::H;
                }
            }
        }
    }

    AlignmentResult {
        score: best_score,
        aligned_len,
        matches,
    }
}

enum State {
    H,
    E,
    F,
}

fn zero() -> AlignmentResult {
    AlignmentResult {
        score: 0,
        aligned_len: 0,
        matches: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::smith_waterman_gotoh_banded_tb;
    use crate::alignment::gapped_scalar::smith_waterman_gotoh_scalar;
    use crate::alignment::scoring::ScoringParams;

    const AA: &[u8] = b"ARNDCQEGHILKMFPSTWYV";

    fn rand_seq(seed: u64, len: usize) -> Vec<u8> {
        let mut x = seed.wrapping_mul(0x9e3779b97f4a7c15).wrapping_add(0x1234_5678);
        (0..len)
            .map(|_| {
                x ^= x >> 12;
                x ^= x << 25;
                x ^= x >> 27;
                AA[(x >> 33) as usize % AA.len()]
            })
            .collect()
    }

    fn mutate(base: &[u8], seed: u64, sub_pct: u64) -> Vec<u8> {
        let mut x = seed | 1;
        base.iter()
            .map(|&c| {
                x ^= x >> 12;
                x ^= x << 25;
                x ^= x >> 27;
                if x % 100 < sub_pct {
                    AA[(x >> 20) as usize % AA.len()]
                } else {
                    c
                }
            })
            .collect()
    }

    #[test]
    fn traceback_score_matches_reference() {
        let scorings = [
            ScoringParams::protein_default(),
            ScoringParams::protein_blosum62(),
        ];
        for s in &scorings {
            for seed in 0..40u64 {
                let la = 20 + (seed as usize % 140);
                let a = rand_seq(seed, la);
                let mutated = mutate(&a, seed.wrapping_mul(11).wrapping_add(5), 30);
                let n = 8 + (seed as usize % la.max(1));
                let b = &mutated[..n.min(mutated.len())];

                let reference = smith_waterman_gotoh_scalar(&a, b, s);
                let tb = smith_waterman_gotoh_banded_tb(&a, b, s, a.len() + b.len());
                // Full band => optimal score must match the reference exactly.
                assert_eq!(tb.score, reference.score, "score seed={seed}");
                // Traceback identity must be self-consistent.
                assert!(tb.matches <= tb.aligned_len, "tb={tb:?}");
                if reference.score > 0 {
                    assert!(tb.aligned_len > 0);
                }
            }
        }
    }

    #[test]
    fn identical_sequences_are_fully_aligned() {
        let s = ScoringParams::protein_blosum62();
        let a = rand_seq(7, 100);
        let r = smith_waterman_gotoh_banded_tb(&a, &a, &s, 200);
        assert_eq!(r.matches, 100);
        assert_eq!(r.aligned_len, 100);
    }

    #[test]
    fn small_offset_alignment_is_captured() {
        let s = ScoringParams::protein_blosum62();
        let a = rand_seq(99, 80);
        let mut b = vec![b'Q', b'Q'];
        b.extend_from_slice(&a);
        let r = smith_waterman_gotoh_banded_tb(&a, &b, &s, 16);
        assert!(r.matches >= 78, "{r:?}");
    }
}
