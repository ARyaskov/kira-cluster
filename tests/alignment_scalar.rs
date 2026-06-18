#[path = "../src/alignment/mod.rs"]
mod alignment;
#[path = "../src/error.rs"]
mod error;
#[path = "../src/simd/mod.rs"]
mod simd;

use alignment::gapped_scalar::smith_waterman_gotoh_scalar;
use alignment::scoring::ScoringParams;

#[test]
fn scalar_alignment_finds_local_match() {
    let a = b"AAAAACCCCCGGGGG";
    let b = b"TTTTACCCCCGGGGGTTTT";
    let scoring = ScoringParams::protein_default();

    let r = smith_waterman_gotoh_scalar(a, b, &scoring);
    assert!(r.score > 0);
    assert!(r.matches >= 10, "{r:?}");
    assert_eq!(r.aligned_len, r.matches, "{r:?}");
}
