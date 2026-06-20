use crate::alignment::AlignmentResult;
#[cfg(feature = "verify-simd")]
use crate::alignment::gapped_scalar::smith_waterman_gotoh_scalar;
use crate::alignment::scoring::ScoringParams;

const NEG_INF: i32 = i32::MIN / 4;
const MAX_LANES: usize = 8;

#[derive(Clone, Copy)]
struct Cell {
    score: i32,
    aligned_len: u32,
    matches: u32,
}

impl Cell {
    const ZERO: Self = Self {
        score: 0,
        aligned_len: 0,
        matches: 0,
    };

    const NEG: Self = Self {
        score: NEG_INF,
        aligned_len: 0,
        matches: 0,
    };
}

#[derive(Clone, Copy)]
enum ScoreKernelBackend {
    Avx2,
    Neon,
}

struct Diagonal {
    sum: usize,
    start_i: usize,
    h: Vec<Cell>,
    e: Vec<Cell>,
    f: Vec<Cell>,
}

impl Diagonal {
    fn empty(sum: usize) -> Self {
        Self {
            sum,
            start_i: 1,
            h: Vec::new(),
            e: Vec::new(),
            f: Vec::new(),
        }
    }

    fn new(sum: usize, start_i: usize, len: usize) -> Self {
        Self {
            sum,
            start_i,
            h: vec![Cell::ZERO; len],
            e: vec![Cell::NEG; len],
            f: vec![Cell::NEG; len],
        }
    }

    fn cell_at(cells: &[Cell], sum: usize, start_i: usize, i: usize, j: usize) -> Cell {
        if i == 0 || j == 0 {
            return Cell::ZERO;
        }
        if i + j != sum || i < start_i {
            return Cell::NEG;
        }
        cells.get(i - start_i).copied().unwrap_or(Cell::NEG)
    }

    fn gap_at(cells: &[Cell], sum: usize, start_i: usize, i: usize, j: usize) -> Cell {
        if i == 0 || j == 0 {
            return Cell::NEG;
        }
        Self::cell_at(cells, sum, start_i, i, j)
    }

    fn h_at(&self, i: usize, j: usize) -> Cell {
        Self::cell_at(&self.h, self.sum, self.start_i, i, j)
    }

    fn e_at(&self, i: usize, j: usize) -> Cell {
        Self::gap_at(&self.e, self.sum, self.start_i, i, j)
    }

    fn f_at(&self, i: usize, j: usize) -> Cell {
        Self::gap_at(&self.f, self.sum, self.start_i, i, j)
    }
}

#[derive(Clone, Copy)]
struct ScoreInputs {
    a: [u8; MAX_LANES],
    b: [u8; MAX_LANES],
    diag: [i32; MAX_LANES],
    h_left: [i32; MAX_LANES],
    e_left: [i32; MAX_LANES],
    h_up: [i32; MAX_LANES],
    f_up: [i32; MAX_LANES],
}

impl ScoreInputs {
    fn new() -> Self {
        Self {
            a: [0; MAX_LANES],
            b: [0; MAX_LANES],
            diag: [NEG_INF; MAX_LANES],
            h_left: [NEG_INF; MAX_LANES],
            e_left: [NEG_INF; MAX_LANES],
            h_up: [NEG_INF; MAX_LANES],
            f_up: [NEG_INF; MAX_LANES],
        }
    }
}

#[derive(Clone, Copy)]
struct ScoreOutputs {
    sub: [i32; MAX_LANES],
    e: [i32; MAX_LANES],
    f: [i32; MAX_LANES],
    h: [i32; MAX_LANES],
}

impl ScoreOutputs {
    fn new() -> Self {
        Self {
            sub: [0; MAX_LANES],
            e: [NEG_INF; MAX_LANES],
            f: [NEG_INF; MAX_LANES],
            h: [0; MAX_LANES],
        }
    }
}

pub fn smith_waterman_gotoh_simd_avx2(
    a: &[u8],
    b: &[u8],
    scoring: &ScoringParams,
) -> AlignmentResult {
    smith_waterman_gotoh_antidiagonal(a, b, scoring, ScoreKernelBackend::Avx2)
}

pub fn smith_waterman_gotoh_simd_neon(
    a: &[u8],
    b: &[u8],
    scoring: &ScoringParams,
) -> AlignmentResult {
    smith_waterman_gotoh_antidiagonal(a, b, scoring, ScoreKernelBackend::Neon)
}

fn smith_waterman_gotoh_antidiagonal(
    a: &[u8],
    b: &[u8],
    scoring: &ScoringParams,
    backend: ScoreKernelBackend,
) -> AlignmentResult {
    if a.is_empty() || b.is_empty() {
        return AlignmentResult {
            score: 0,
            aligned_len: 0,
            matches: 0,
        };
    }

    let m = a.len();
    let n = b.len();
    let go_ge = scoring.gap_open as i32 + scoring.gap_extend as i32;
    let ge = scoring.gap_extend as i32;
    let lanes = lanes_for_backend(backend);

    let mut prev2 = Diagonal::empty(0);
    let mut prev1 = Diagonal::empty(1);
    let mut best = Cell::ZERO;

    for sum in 2..=m + n {
        let Some((start_i, end_i)) = diagonal_bounds(sum, m, n) else {
            prev2 = prev1;
            prev1 = Diagonal::empty(sum);
            continue;
        };
        let len = end_i - start_i + 1;
        let mut curr = Diagonal::new(sum, start_i, len);

        let mut base = 0usize;
        while base < len {
            let chunk_len = (len - base).min(lanes);
            let mut inputs = ScoreInputs::new();
            let mut diag_cells = [Cell::NEG; MAX_LANES];
            let mut left_h_cells = [Cell::NEG; MAX_LANES];
            let mut left_e_cells = [Cell::NEG; MAX_LANES];
            let mut up_h_cells = [Cell::NEG; MAX_LANES];
            let mut up_f_cells = [Cell::NEG; MAX_LANES];

            for lane in 0..chunk_len {
                let i = start_i + base + lane;
                let j = sum - i;
                let diag = prev2.h_at(i - 1, j - 1);
                let left_h = prev1.h_at(i, j - 1);
                let left_e = prev1.e_at(i, j - 1);
                let up_h = prev1.h_at(i - 1, j);
                let up_f = prev1.f_at(i - 1, j);

                diag_cells[lane] = diag;
                left_h_cells[lane] = left_h;
                left_e_cells[lane] = left_e;
                up_h_cells[lane] = up_h;
                up_f_cells[lane] = up_f;

                inputs.a[lane] = a[i - 1];
                inputs.b[lane] = b[j - 1];
                inputs.diag[lane] = diag.score;
                inputs.h_left[lane] = left_h.score;
                inputs.e_left[lane] = left_e.score;
                inputs.h_up[lane] = up_h.score;
                inputs.f_up[lane] = up_f.score;
            }

            let scores = score_chunk(backend, chunk_len, &inputs, scoring, go_ge, ge);

            for lane in 0..chunk_len {
                let idx = base + lane;
                let is_match = inputs.a[lane] == inputs.b[lane];

                let diag = add_score(diag_cells[lane], scores.sub[lane], is_match);

                let mut e = pick_better(
                    add_gap(left_h_cells[lane], go_ge),
                    add_gap(left_e_cells[lane], ge),
                );
                e.score = normalize_gap_score(scores.e[lane]);

                let mut f = pick_better(
                    add_gap(up_h_cells[lane], go_ge),
                    add_gap(up_f_cells[lane], ge),
                );
                f.score = normalize_gap_score(scores.f[lane]);

                let mut h = pick_best_local(diag, e, f);
                if scores.h[lane] <= 0 {
                    h = Cell::ZERO;
                } else {
                    h.score = scores.h[lane];
                }

                curr.e[idx] = e;
                curr.f[idx] = f;
                curr.h[idx] = h;
                best = pick_better(best, h);
            }

            base += chunk_len;
        }

        prev2 = prev1;
        prev1 = curr;
    }

    let result = AlignmentResult {
        score: best.score,
        aligned_len: best.aligned_len,
        matches: best.matches,
    };

    #[cfg(feature = "verify-simd")]
    debug_assert_eq!(result, smith_waterman_gotoh_scalar(a, b, scoring));
    result
}

fn diagonal_bounds(sum: usize, m: usize, n: usize) -> Option<(usize, usize)> {
    let start = sum.saturating_sub(n).max(1);
    let end = (sum - 1).min(m);
    if start <= end {
        Some((start, end))
    } else {
        None
    }
}

fn lanes_for_backend(backend: ScoreKernelBackend) -> usize {
    match backend {
        ScoreKernelBackend::Avx2 => 8,
        ScoreKernelBackend::Neon => 4,
    }
}

fn normalize_gap_score(score: i32) -> i32 {
    if score <= NEG_INF / 2 { NEG_INF } else { score }
}

fn score_chunk(
    backend: ScoreKernelBackend,
    lanes: usize,
    inputs: &ScoreInputs,
    scoring: &ScoringParams,
    go_ge: i32,
    ge: i32,
) -> ScoreOutputs {
    match backend {
        ScoreKernelBackend::Avx2 => score_chunk_avx2(lanes, inputs, scoring, go_ge, ge),
        ScoreKernelBackend::Neon => score_chunk_neon(lanes, inputs, scoring, go_ge, ge),
    }
}

fn fill_sub_scores(
    lanes: usize,
    inputs: &ScoreInputs,
    scoring: &ScoringParams,
    out: &mut ScoreOutputs,
) {
    // Substitution scores are computed scalar per lane (cheap table lookup or
    // match/mismatch); the surrounding H/E/F recurrence is what gets vectorized.
    // Routing through `score_pair` keeps the SIMD path identical to the scalar
    // path for both flat and BLOSUM62 scoring.
    for lane in 0..lanes {
        out.sub[lane] = scoring.score_pair(inputs.a[lane], inputs.b[lane]) as i32;
    }
}

fn score_chunk_scalar(
    lanes: usize,
    inputs: &ScoreInputs,
    scoring: &ScoringParams,
    go_ge: i32,
    ge: i32,
) -> ScoreOutputs {
    let mut out = ScoreOutputs::new();
    fill_sub_scores(lanes, inputs, scoring, &mut out);
    for lane in 0..lanes {
        let diag = inputs.diag[lane] + out.sub[lane];
        let e = (inputs.h_left[lane] + go_ge).max(inputs.e_left[lane] + ge);
        let f = (inputs.h_up[lane] + go_ge).max(inputs.f_up[lane] + ge);
        out.e[lane] = e;
        out.f[lane] = f;
        out.h[lane] = 0.max(diag).max(e).max(f);
    }
    out
}

#[cfg(target_arch = "x86_64")]
fn score_chunk_avx2(
    lanes: usize,
    inputs: &ScoreInputs,
    scoring: &ScoringParams,
    go_ge: i32,
    ge: i32,
) -> ScoreOutputs {
    if std::arch::is_x86_feature_detected!("avx2") {
        // SAFETY: guarded by runtime AVX2 detection; all loads/stores use fixed
        // stack arrays with at least eight i32 lanes.
        unsafe { score_chunk_avx2_impl(lanes, inputs, scoring, go_ge, ge) }
    } else {
        score_chunk_scalar(lanes, inputs, scoring, go_ge, ge)
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn score_chunk_avx2_impl(
    lanes: usize,
    inputs: &ScoreInputs,
    scoring: &ScoringParams,
    go_ge: i32,
    ge: i32,
) -> ScoreOutputs {
    use core::arch::x86_64::{
        __m256i, _mm256_add_epi32, _mm256_loadu_si256, _mm256_max_epi32, _mm256_set1_epi32,
        _mm256_storeu_si256,
    };

    let mut out = ScoreOutputs::new();
    fill_sub_scores(lanes, inputs, scoring, &mut out);

    let diag_v = unsafe { _mm256_loadu_si256(inputs.diag.as_ptr() as *const __m256i) };
    let sub_v = unsafe { _mm256_loadu_si256(out.sub.as_ptr() as *const __m256i) };
    let h_left_v = unsafe { _mm256_loadu_si256(inputs.h_left.as_ptr() as *const __m256i) };
    let e_left_v = unsafe { _mm256_loadu_si256(inputs.e_left.as_ptr() as *const __m256i) };
    let h_up_v = unsafe { _mm256_loadu_si256(inputs.h_up.as_ptr() as *const __m256i) };
    let f_up_v = unsafe { _mm256_loadu_si256(inputs.f_up.as_ptr() as *const __m256i) };

    let go_ge_v = _mm256_set1_epi32(go_ge);
    let ge_v = _mm256_set1_epi32(ge);
    let zero_v = _mm256_set1_epi32(0);

    let diag_score = _mm256_add_epi32(diag_v, sub_v);
    let e_score = _mm256_max_epi32(
        _mm256_add_epi32(h_left_v, go_ge_v),
        _mm256_add_epi32(e_left_v, ge_v),
    );
    let f_score = _mm256_max_epi32(
        _mm256_add_epi32(h_up_v, go_ge_v),
        _mm256_add_epi32(f_up_v, ge_v),
    );
    let h_score = _mm256_max_epi32(
        zero_v,
        _mm256_max_epi32(diag_score, _mm256_max_epi32(e_score, f_score)),
    );

    unsafe { _mm256_storeu_si256(out.e.as_mut_ptr() as *mut __m256i, e_score) };
    unsafe { _mm256_storeu_si256(out.f.as_mut_ptr() as *mut __m256i, f_score) };
    unsafe { _mm256_storeu_si256(out.h.as_mut_ptr() as *mut __m256i, h_score) };
    out
}

#[cfg(not(target_arch = "x86_64"))]
fn score_chunk_avx2(
    lanes: usize,
    inputs: &ScoreInputs,
    scoring: &ScoringParams,
    go_ge: i32,
    ge: i32,
) -> ScoreOutputs {
    score_chunk_scalar(lanes, inputs, scoring, go_ge, ge)
}

#[cfg(target_arch = "aarch64")]
fn score_chunk_neon(
    lanes: usize,
    inputs: &ScoreInputs,
    scoring: &ScoringParams,
    go_ge: i32,
    ge: i32,
) -> ScoreOutputs {
    let mut out = ScoreOutputs::new();
    let mut base = 0usize;
    while base < lanes {
        let width = (lanes - base).min(4);
        // SAFETY: NEON is mandatory on aarch64 and the helper loads from fixed
        // stack arrays with at least four i32 lanes from the requested base.
        unsafe { score_chunk_neon_impl(base, width, inputs, scoring, go_ge, ge, &mut out) };
        base += width;
    }
    out
}

#[cfg(target_arch = "aarch64")]
unsafe fn score_chunk_neon_impl(
    base: usize,
    lanes: usize,
    inputs: &ScoreInputs,
    scoring: &ScoringParams,
    go_ge: i32,
    ge: i32,
    out: &mut ScoreOutputs,
) {
    use core::arch::aarch64::{vaddq_s32, vdupq_n_s32, vld1q_s32, vmaxq_s32, vst1q_s32};

    fill_sub_scores(base + lanes, inputs, scoring, out);

    let diag_v = unsafe { vld1q_s32(inputs.diag[base..].as_ptr()) };
    let sub_v = unsafe { vld1q_s32(out.sub[base..].as_ptr()) };
    let h_left_v = unsafe { vld1q_s32(inputs.h_left[base..].as_ptr()) };
    let e_left_v = unsafe { vld1q_s32(inputs.e_left[base..].as_ptr()) };
    let h_up_v = unsafe { vld1q_s32(inputs.h_up[base..].as_ptr()) };
    let f_up_v = unsafe { vld1q_s32(inputs.f_up[base..].as_ptr()) };

    let go_ge_v = unsafe { vdupq_n_s32(go_ge) };
    let ge_v = unsafe { vdupq_n_s32(ge) };
    let zero_v = unsafe { vdupq_n_s32(0) };

    let diag_score = unsafe { vaddq_s32(diag_v, sub_v) };
    let e_score = unsafe { vmaxq_s32(vaddq_s32(h_left_v, go_ge_v), vaddq_s32(e_left_v, ge_v)) };
    let f_score = unsafe { vmaxq_s32(vaddq_s32(h_up_v, go_ge_v), vaddq_s32(f_up_v, ge_v)) };
    let h_score = unsafe { vmaxq_s32(zero_v, vmaxq_s32(diag_score, vmaxq_s32(e_score, f_score))) };

    unsafe { vst1q_s32(out.e[base..].as_mut_ptr(), e_score) };
    unsafe { vst1q_s32(out.f[base..].as_mut_ptr(), f_score) };
    unsafe { vst1q_s32(out.h[base..].as_mut_ptr(), h_score) };
}

#[cfg(not(target_arch = "aarch64"))]
fn score_chunk_neon(
    lanes: usize,
    inputs: &ScoreInputs,
    scoring: &ScoringParams,
    go_ge: i32,
    ge: i32,
) -> ScoreOutputs {
    score_chunk_scalar(lanes, inputs, scoring, go_ge, ge)
}

fn add_score(src: Cell, score: i32, is_match: bool) -> Cell {
    if src.score <= NEG_INF / 2 {
        return Cell::NEG;
    }
    Cell {
        score: src.score + score,
        aligned_len: src.aligned_len + 1,
        matches: src.matches + u32::from(is_match),
    }
}

fn add_gap(src: Cell, gap_score: i32) -> Cell {
    if src.score <= NEG_INF / 2 {
        return Cell::NEG;
    }
    Cell {
        score: src.score + gap_score,
        aligned_len: src.aligned_len + 1,
        matches: src.matches,
    }
}

fn pick_best_local(a: Cell, b: Cell, c: Cell) -> Cell {
    let best_nonzero = pick_better(pick_better(a, b), c);
    if best_nonzero.score <= 0 {
        Cell::ZERO
    } else {
        best_nonzero
    }
}

fn pick_better(a: Cell, b: Cell) -> Cell {
    if b.score > a.score {
        return b;
    }
    if b.score < a.score {
        return a;
    }
    if b.matches > a.matches {
        return b;
    }
    if b.matches < a.matches {
        return a;
    }
    if b.aligned_len < a.aligned_len {
        return b;
    }
    a
}
