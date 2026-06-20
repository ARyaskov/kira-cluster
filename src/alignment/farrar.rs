//! Striped (Farrar) SIMD Smith-Waterman **score** kernel.
//!
//! Computes only the optimal local-alignment score — no traceback — which is the
//! part that vectorizes well. It is used as a fast pre-screen: pairs whose score
//! is far too low to possibly meet the identity/coverage gate skip the expensive
//! traceback DP entirely (MMseqs2-style ungapped/gapped prefiltering).
//!
//! The AVX2 path processes 16 query residues per vector with a precomputed query
//! profile (so there is no per-cell substitution-matrix lookup) and the standard
//! lazy-F correction. i16 saturation is detected and falls back to the scalar
//! score so the result is never silently wrong.

use crate::alignment::scoring::ScoringParams;
use crate::simd::SimdBackend;

/// Score-only local Gotoh alignment. Returns the optimal SW score (>= 0).
pub fn sw_score(a: &[u8], b: &[u8], scoring: &ScoringParams, backend: SimdBackend) -> i32 {
    if a.is_empty() || b.is_empty() {
        return 0;
    }
    match backend {
        SimdBackend::Avx2 => sw_score_avx2(a, b, scoring),
        _ => sw_score_scalar(a, b, scoring),
    }
}

/// Reference score-only DP (also the saturation fallback). O(|a|*|b|), O(|b|)
/// memory. Identical optimal score to `gapped_scalar::smith_waterman_gotoh_scalar`.
pub fn sw_score_scalar(a: &[u8], b: &[u8], scoring: &ScoringParams) -> i32 {
    let n = b.len();
    let go_ge = scoring.gap_open as i32 + scoring.gap_extend as i32; // negative
    let ge = scoring.gap_extend as i32; // negative
    let neg = i32::MIN / 4;

    let mut prev_h = vec![0i32; n + 1];
    let mut curr_h = vec![0i32; n + 1];
    let mut prev_f = vec![neg; n + 1];
    let mut curr_f = vec![neg; n + 1];
    let mut best = 0i32;

    for i in 1..=a.len() {
        let ai = a[i - 1];
        let mut e = neg;
        curr_h[0] = 0;
        curr_f[0] = neg;
        for j in 1..=n {
            let sub = scoring.score_pair(ai, b[j - 1]) as i32;
            let diag = prev_h[j - 1] + sub;
            e = (curr_h[j - 1] + go_ge).max(e + ge);
            curr_f[j] = (prev_h[j] + go_ge).max(prev_f[j] + ge);
            let h = 0.max(diag).max(e).max(curr_f[j]);
            curr_h[j] = h;
            if h > best {
                best = h;
            }
        }
        std::mem::swap(&mut prev_h, &mut curr_h);
        std::mem::swap(&mut prev_f, &mut curr_f);
    }
    best
}

#[cfg(target_arch = "x86_64")]
fn sw_score_avx2(a: &[u8], b: &[u8], scoring: &ScoringParams) -> i32 {
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: guarded by runtime AVX2 detection.
        match unsafe { sw_score_avx2_impl(a, b, scoring) } {
            Some(score) => score,
            None => sw_score_scalar(a, b, scoring), // i16 saturated
        }
    } else {
        sw_score_scalar(a, b, scoring)
    }
}

#[cfg(not(target_arch = "x86_64"))]
fn sw_score_avx2(a: &[u8], b: &[u8], scoring: &ScoringParams) -> i32 {
    sw_score_scalar(a, b, scoring)
}

#[cfg(target_arch = "x86_64")]
const LANES: usize = 16; // i16 lanes in a 256-bit register

/// Striped Farrar kernel. Returns `None` if the score saturated i16 (caller
/// falls back to the scalar score). Query = `a`, target = `b`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn sw_score_avx2_impl(a: &[u8], b: &[u8], scoring: &ScoringParams) -> Option<i32> {
    use core::arch::x86_64::*;

    let qlen = a.len();
    let seg_len = qlen.div_ceil(LANES); // segments per column
    if seg_len == 0 {
        return Some(0);
    }

    // Gap penalties as positive magnitudes (scores are subtracted).
    let gap_open_ext = -(scoring.gap_open as i32 + scoring.gap_extend as i32); // first gapped cell
    let gap_ext = -(scoring.gap_extend as i32);
    if gap_open_ext > i16::MAX as i32 || gap_ext > i16::MAX as i32 {
        return None;
    }
    let v_gap_oe = _mm256_set1_epi16(gap_open_ext as i16);
    let v_gap_e = _mm256_set1_epi16(gap_ext as i16);
    let v_zero = _mm256_setzero_si256();

    // Build the striped query profile, indexed by BLOSUM residue class.
    // profile[class][seg] : lane k holds score(query[k*seg_len + seg], class).
    let num_classes = 24usize;
    let mut profile: Vec<__m256i> = vec![v_zero; num_classes * seg_len];
    {
        let mut tmp = [0i16; LANES];
        for class in 0..num_classes {
            let class_byte = class_to_byte(class);
            for seg in 0..seg_len {
                for (k, slot) in tmp.iter_mut().enumerate() {
                    let qi = k * seg_len + seg;
                    *slot = if qi < qlen {
                        scoring.score_pair(a[qi], class_byte)
                    } else {
                        // Padding query positions: large negative so they never win.
                        i16::MIN / 2
                    };
                }
                // SAFETY: tmp has exactly LANES i16 elements.
                profile[class * seg_len + seg] =
                    unsafe { _mm256_loadu_si256(tmp.as_ptr() as *const __m256i) };
            }
        }
    }

    let mut pv_h_store: Vec<__m256i> = vec![v_zero; seg_len];
    let mut pv_h_load: Vec<__m256i> = vec![v_zero; seg_len];
    let mut pv_e: Vec<__m256i> = vec![_mm256_set1_epi16(i16::MIN / 2); seg_len];

    let sat = _mm256_set1_epi16(i16::MAX);
    let mut v_max = v_zero;

    for &bj in b {
        let class = byte_to_class(bj);
        let prof = &profile[class * seg_len..(class + 1) * seg_len];

        // vH initialized from the last segment of the previous column, shifted
        // one i16 lane towards higher query positions (cross-segment carry).
        // SAFETY: avx2 is enabled for this function.
        let mut v_h = unsafe { shift_left_epi16(pv_h_store[seg_len - 1]) };
        let mut v_f = v_zero;

        std::mem::swap(&mut pv_h_load, &mut pv_h_store);

        for seg in 0..seg_len {
            v_h = _mm256_adds_epi16(v_h, prof[seg]);
            v_h = _mm256_max_epi16(v_h, pv_e[seg]);
            v_h = _mm256_max_epi16(v_h, v_f);
            v_h = _mm256_max_epi16(v_h, v_zero);
            v_max = _mm256_max_epi16(v_max, v_h);
            pv_h_store[seg] = v_h;

            // E and F for the next cell.
            let h_minus_open = _mm256_subs_epi16(v_h, v_gap_oe);
            pv_e[seg] = _mm256_max_epi16(_mm256_subs_epi16(pv_e[seg], v_gap_e), h_minus_open);
            v_f = _mm256_max_epi16(_mm256_subs_epi16(v_f, v_gap_e), h_minus_open);

            v_h = pv_h_load[seg];
        }

        // Lazy-F loop: propagate F across segment boundaries until no lane of F
        // can open a new gap (Farrar's correction).
        'lazy: for _ in 0..LANES {
            // SAFETY: avx2 is enabled for this function.
            v_f = unsafe { shift_left_epi16(v_f) };
            for h_slot in pv_h_store.iter_mut() {
                let h_new = _mm256_max_epi16(*h_slot, v_f);
                *h_slot = h_new;
                v_max = _mm256_max_epi16(v_max, h_new);
                v_f = _mm256_subs_epi16(v_f, v_gap_e);
                // Done once no lane has F > H - gap_open_ext (F can't open a gap).
                let thresh = _mm256_subs_epi16(h_new, v_gap_oe);
                let gt = _mm256_cmpgt_epi16(v_f, thresh);
                if _mm256_testz_si256(gt, gt) == 1 {
                    break 'lazy;
                }
            }
        }
    }

    // Horizontal max over v_max (16 i16 lanes).
    // SAFETY: avx2 is enabled for this function.
    let best = unsafe { hmax_epi16(v_max) } as i32;
    // Saturation check: if we hit i16::MAX anywhere the result is unreliable.
    let hit_sat = _mm256_cmpeq_epi16(v_max, sat);
    if _mm256_testz_si256(hit_sat, hit_sat) == 0 {
        return None;
    }
    Some(best)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn shift_left_epi16(v: core::arch::x86_64::__m256i) -> core::arch::x86_64::__m256i {
    use core::arch::x86_64::*;
    // Shift the 256-bit value left by one i16 lane across the full register.
    // _mm256_slli_si256 operates per 128-bit lane, so we also carry the top
    // element of the low lane into the high lane.
    let shifted = _mm256_slli_si256(v, 2);
    // Bring element 7 (top of low 128) into element 8 (bottom of high 128).
    let carry = _mm256_permute2x128_si256(v, v, 0x08); // low lane -> high position, zero low
    let carry = _mm256_srli_si256(carry, 14); // top i16 of low lane to lane-8 slot
    _mm256_or_si256(shifted, carry)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn hmax_epi16(v: core::arch::x86_64::__m256i) -> i16 {
    use core::arch::x86_64::*;
    let mut buf = [0i16; LANES];
    // SAFETY: buf has exactly LANES i16 elements.
    unsafe { _mm256_storeu_si256(buf.as_mut_ptr() as *mut __m256i, v) };
    let mut m = i16::MIN;
    for x in buf {
        if x > m {
            m = x;
        }
    }
    m
}

/// BLOSUM residue-class index for a byte (mirrors scoring::BLOSUM_INDEX order).
fn byte_to_class(b: u8) -> usize {
    const ORDER: &[u8; 23] = b"ARNDCQEGHILKMFPSTWYVBZX";
    let up = b.to_ascii_uppercase();
    ORDER.iter().position(|&c| c == up).unwrap_or(22)
}

fn class_to_byte(class: usize) -> u8 {
    const ORDER: &[u8; 23] = b"ARNDCQEGHILKMFPSTWYVBZX";
    *ORDER.get(class).unwrap_or(&b'X')
}

#[cfg(test)]
mod tests {
    use super::{sw_score_avx2, sw_score_scalar};
    use crate::alignment::gapped_scalar::smith_waterman_gotoh_scalar;
    use crate::alignment::scoring::ScoringParams;

    const AA: &[u8] = b"ARNDCQEGHILKMFPSTWYV";

    fn rand_seq(seed: u64, len: usize) -> Vec<u8> {
        let mut x = seed.wrapping_mul(0x9e3779b97f4a7c15).wrapping_add(0x55);
        (0..len)
            .map(|_| {
                x ^= x >> 12;
                x ^= x << 25;
                x ^= x >> 27;
                AA[(x >> 33) as usize % AA.len()]
            })
            .collect()
    }

    #[test]
    fn scalar_score_matches_reference() {
        let s = ScoringParams::protein_blosum62();
        for seed in 0..40u64 {
            let a = rand_seq(seed, 15 + (seed as usize % 90));
            let b = rand_seq(seed.wrapping_mul(3) + 1, 15 + (seed as usize % 110));
            let want = smith_waterman_gotoh_scalar(&a, &b, &s).score;
            assert_eq!(sw_score_scalar(&a, &b, &s), want, "seed={seed}");
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_score_matches_reference() {
        let scorings = [
            ScoringParams::protein_blosum62(),
            ScoringParams::protein_default(),
        ];
        for s in &scorings {
            for seed in 0..80u64 {
                let a = rand_seq(seed, 10 + (seed as usize % 130));
                let b = rand_seq(seed.wrapping_mul(7) + 5, 10 + (seed as usize % 150));
                let want = smith_waterman_gotoh_scalar(&a, &b, s).score;
                let got = sw_score_avx2(&a, &b, s);
                assert_eq!(got, want, "seed={seed}, a.len={}, b.len={}", a.len(), b.len());
            }
        }
    }
}
