use crate::simd::SimdBackend;

pub fn hamming_filter(a: &[u8], b: &[u8], max_mismatches: u32, backend: SimdBackend) -> bool {
    let len = a.len().min(b.len());
    let a = &a[..len];
    let b = &b[..len];

    let mismatches = match backend {
        SimdBackend::Avx2 => hamming_count_avx2(a, b, max_mismatches),
        SimdBackend::Neon => hamming_count_neon(a, b, max_mismatches),
        SimdBackend::Scalar => hamming_count_scalar(a, b, max_mismatches),
    };

    mismatches <= max_mismatches
}

pub fn ungapped_filter(
    a: &[u8],
    b: &[u8],
    min_identity: f32,
    min_coverage: f32,
    backend: SimdBackend,
) -> bool {
    let overlap = a.len().min(b.len());
    if overlap == 0 {
        return false;
    }

    let coverage = (overlap as f32) / (a.len().max(b.len()) as f32);
    if coverage < min_coverage {
        return false;
    }

    let allowed = ((overlap as f32) * (1.0 - min_identity)).floor().max(0.0) as u32;
    let mismatches = match backend {
        SimdBackend::Avx2 => hamming_count_avx2(&a[..overlap], &b[..overlap], u32::MAX),
        SimdBackend::Neon => hamming_count_neon(&a[..overlap], &b[..overlap], u32::MAX),
        SimdBackend::Scalar => hamming_count_scalar(&a[..overlap], &b[..overlap], u32::MAX),
    };

    if mismatches > allowed {
        return false;
    }

    let matches = overlap as u32 - mismatches;
    let identity = (matches as f32) / (overlap as f32);
    identity >= min_identity
}

/// MMseqs2-style coverage decision. `q` is the query (the sequence being
/// assigned / searched), `t` is the target (the cluster representative / DB hit).
///
/// - mode 0: bidirectional — alignment must cover `min_cov` of BOTH query and
///   target (equivalently `aligned_len / max(q, t) >= min_cov`).
/// - mode 1: coverage of target (`aligned_len / t`).
/// - mode 2: coverage of query (`aligned_len / q`).
/// - mode 3: target length is at least `min_cov` of the query length (a pure
///   length-ratio gate; the alignment length is ignored).
pub fn coverage_satisfied(
    aligned_len: u32,
    qlen: usize,
    tlen: usize,
    min_cov: f32,
    cov_mode: u8,
) -> bool {
    let a = aligned_len as f32;
    let q = qlen.max(1) as f32;
    let t = tlen.max(1) as f32;
    match cov_mode {
        1 => a / t >= min_cov,
        2 => a / q >= min_cov,
        3 => t / q >= min_cov,
        _ => a / q >= min_cov && a / t >= min_cov,
    }
}

/// Length-only necessary condition for [`coverage_satisfied`]: the best possible
/// coverage is achieved when the alignment spans the entire shorter sequence
/// (`aligned_len <= min(q, t)`). If even that cannot reach `min_cov`, the pair
/// can be rejected before alignment without losing any true positive.
pub fn coverage_prefilter(qlen: usize, tlen: usize, min_cov: f32, cov_mode: u8) -> bool {
    let best = qlen.min(tlen) as u32;
    coverage_satisfied(best, qlen, tlen, min_cov, cov_mode)
}

fn hamming_count_scalar(a: &[u8], b: &[u8], max_mismatches: u32) -> u32 {
    let mut mismatches = 0u32;
    for (&x, &y) in a.iter().zip(b) {
        if x != y {
            mismatches += 1;
            if mismatches > max_mismatches {
                return mismatches;
            }
        }
    }
    mismatches
}

#[cfg(target_arch = "x86_64")]
fn hamming_count_avx2(a: &[u8], b: &[u8], max_mismatches: u32) -> u32 {
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: guarded by runtime feature detection.
        return unsafe { hamming_count_avx2_impl(a, b, max_mismatches) };
    }
    hamming_count_scalar(a, b, max_mismatches)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn hamming_count_avx2_impl(a: &[u8], b: &[u8], max_mismatches: u32) -> u32 {
    use core::arch::x86_64::{
        __m256i, _mm256_cmpeq_epi8, _mm256_loadu_si256, _mm256_movemask_epi8,
    };

    let mut mismatches = 0u32;
    let mut i = 0usize;
    while i + 32 <= a.len() {
        // SAFETY: bounds checked by loop condition.
        let va = unsafe { _mm256_loadu_si256(a.as_ptr().add(i) as *const __m256i) };
        // SAFETY: bounds checked by loop condition.
        let vb = unsafe { _mm256_loadu_si256(b.as_ptr().add(i) as *const __m256i) };
        let eq = _mm256_cmpeq_epi8(va, vb);
        let mask = _mm256_movemask_epi8(eq) as u32;
        mismatches += 32 - mask.count_ones();
        if mismatches > max_mismatches {
            return mismatches;
        }
        i += 32;
    }

    mismatches + hamming_count_scalar(&a[i..], &b[i..], max_mismatches.saturating_sub(mismatches))
}

#[cfg(not(target_arch = "x86_64"))]
fn hamming_count_avx2(a: &[u8], b: &[u8], max_mismatches: u32) -> u32 {
    hamming_count_scalar(a, b, max_mismatches)
}

#[cfg(target_arch = "aarch64")]
fn hamming_count_neon(a: &[u8], b: &[u8], max_mismatches: u32) -> u32 {
    // SAFETY: NEON is mandatory on aarch64 targets.
    unsafe { hamming_count_neon_impl(a, b, max_mismatches) }
}

#[cfg(target_arch = "aarch64")]
unsafe fn hamming_count_neon_impl(a: &[u8], b: &[u8], max_mismatches: u32) -> u32 {
    use core::arch::aarch64::{vaddvq_u8, vceqq_u8, vcntq_u8, vld1q_u8};

    let mut mismatches = 0u32;
    let mut i = 0usize;
    while i + 16 <= a.len() {
        // SAFETY: bounds checked by loop condition.
        let va = unsafe { vld1q_u8(a.as_ptr().add(i)) };
        // SAFETY: bounds checked by loop condition.
        let vb = unsafe { vld1q_u8(b.as_ptr().add(i)) };
        // SAFETY: NEON operations on registers produced in this function.
        let eq = unsafe { vceqq_u8(va, vb) };
        // SAFETY: NEON operations on registers produced in this function.
        let bits = unsafe { vaddvq_u8(vcntq_u8(eq)) };
        let eq_count = (bits / 8) as u32;
        mismatches += 16 - eq_count;
        if mismatches > max_mismatches {
            return mismatches;
        }
        i += 16;
    }

    mismatches + hamming_count_scalar(&a[i..], &b[i..], max_mismatches.saturating_sub(mismatches))
}

#[cfg(not(target_arch = "aarch64"))]
fn hamming_count_neon(a: &[u8], b: &[u8], max_mismatches: u32) -> u32 {
    hamming_count_scalar(a, b, max_mismatches)
}

#[cfg(test)]
mod tests {
    use super::{hamming_filter, ungapped_filter};
    use crate::simd::SimdBackend;

    #[test]
    fn simd_paths_are_callable() {
        let a = b"AAAAAAAAAAAAAAAAAAAA";
        let b = b"AAAAAAAAAAAATAAAAAAA";
        for backend in [SimdBackend::Scalar, SimdBackend::Avx2, SimdBackend::Neon] {
            assert!(hamming_filter(a, b, 2, backend));
            assert!(ungapped_filter(a, b, 0.9, 0.8, backend));
        }
    }
}
