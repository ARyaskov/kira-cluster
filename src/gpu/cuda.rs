use crate::alignment::scoring::ScoringParams;
use crate::alignment::{AlignmentResult, SeqId};
use crate::cluster::SeqDb;
use crate::error::{AppError, ErrorKind, Result};

#[cfg(feature = "cuda")]
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
struct CudaAlignmentResult {
    score: i32,
    aligned_len: u32,
    matches: u32,
}

#[cfg(feature = "cuda")]
unsafe extern "C" {
    fn kira_cuda_available() -> i32;
    fn kira_cuda_sw_batch(
        a_ptrs: *const *const u8,
        a_lens: *const u32,
        b_ptrs: *const *const u8,
        b_lens: *const u32,
        n_pairs: u32,
        match_score: i16,
        mismatch_score: i16,
        gap_open: i16,
        gap_extend: i16,
        out_results: *mut CudaAlignmentResult,
    ) -> i32;
}

#[cfg(feature = "cuda")]
pub fn is_available() -> bool {
    // SAFETY: pure availability probe with no pointers.
    unsafe { kira_cuda_available() == 1 }
}

#[cfg(not(feature = "cuda"))]
pub fn is_available() -> bool {
    false
}

#[cfg(feature = "cuda")]
pub fn align_pairs_cuda(
    db: &SeqDb,
    pairs: &[(SeqId, SeqId)],
    scoring: &ScoringParams,
) -> Result<Vec<AlignmentResult>> {
    let mut a_ptrs = Vec::with_capacity(pairs.len());
    let mut b_ptrs = Vec::with_capacity(pairs.len());
    let mut a_lens = Vec::with_capacity(pairs.len());
    let mut b_lens = Vec::with_capacity(pairs.len());

    for &(a_id, b_id) in pairs {
        let a = db.seq(a_id);
        let b = db.seq(b_id);
        a_ptrs.push(a.as_ptr());
        b_ptrs.push(b.as_ptr());
        a_lens.push(a.len() as u32);
        b_lens.push(b.len() as u32);
    }

    let mut out = vec![CudaAlignmentResult::default(); pairs.len()];
    // SAFETY: all pointer arrays are valid for the duration of the call and lengths match n_pairs.
    let rc = unsafe {
        kira_cuda_sw_batch(
            a_ptrs.as_ptr(),
            a_lens.as_ptr(),
            b_ptrs.as_ptr(),
            b_lens.as_ptr(),
            pairs.len() as u32,
            scoring.match_score,
            scoring.mismatch_score,
            scoring.gap_open,
            scoring.gap_extend,
            out.as_mut_ptr(),
        )
    };

    if rc != 0 {
        return Err(AppError::new(
            ErrorKind::Unsupported,
            format!("CUDA kernel execution failed: rc={rc}"),
        ));
    }

    Ok(out
        .into_iter()
        .map(|r| AlignmentResult {
            score: r.score,
            aligned_len: r.aligned_len,
            matches: r.matches,
        })
        .collect())
}

#[cfg(not(feature = "cuda"))]
pub fn align_pairs_cuda(
    _db: &SeqDb,
    _pairs: &[(SeqId, SeqId)],
    _scoring: &ScoringParams,
) -> Result<Vec<AlignmentResult>> {
    Err(AppError::new(
        ErrorKind::Unsupported,
        "CUDA feature not enabled at compile time",
    ))
}
