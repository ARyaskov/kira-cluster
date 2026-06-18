use crate::alignment::AlignmentResult;
use crate::alignment::scoring::ScoringParams;

const NEG_INF: i32 = i32::MIN / 4;

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

pub fn smith_waterman_gotoh_scalar(a: &[u8], b: &[u8], scoring: &ScoringParams) -> AlignmentResult {
    if a.is_empty() || b.is_empty() {
        return AlignmentResult {
            score: 0,
            aligned_len: 0,
            matches: 0,
        };
    }

    let n = b.len();
    let go_ge = scoring.gap_open as i32 + scoring.gap_extend as i32;
    let ge = scoring.gap_extend as i32;

    let mut prev_h = vec![Cell::ZERO; n + 1];
    let mut curr_h = vec![Cell::ZERO; n + 1];
    let mut prev_f = vec![Cell::NEG; n + 1];
    let mut curr_f = vec![Cell::NEG; n + 1];

    let mut best = Cell::ZERO;

    for i in 1..=a.len() {
        curr_h[0] = Cell::ZERO;
        curr_f[0] = Cell::NEG;
        let mut e = Cell::NEG;

        for j in 1..=n {
            let sub = scoring.score_pair(a[i - 1], b[j - 1]) as i32;

            let diag = add_score(prev_h[j - 1], sub, a[i - 1] == b[j - 1]);

            let e_open = add_gap(curr_h[j - 1], go_ge);
            let e_ext = add_gap(e, ge);
            e = pick_better(e_open, e_ext);

            let f_open = add_gap(prev_h[j], go_ge);
            let f_ext = add_gap(prev_f[j], ge);
            curr_f[j] = pick_better(f_open, f_ext);

            let h = pick_best_local(diag, e, curr_f[j]);
            curr_h[j] = h;
            best = pick_better(best, h);
        }

        std::mem::swap(&mut prev_h, &mut curr_h);
        std::mem::swap(&mut prev_f, &mut curr_f);
    }

    AlignmentResult {
        score: best.score,
        aligned_len: best.aligned_len,
        matches: best.matches,
    }
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
