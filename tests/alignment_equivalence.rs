#[path = "../src/alignment/mod.rs"]
mod alignment;
#[path = "../src/error.rs"]
mod error;
#[path = "../src/simd/mod.rs"]
mod simd;

use alignment::scoring::ScoringParams;
use alignment::{AlignmentResult, smith_waterman_gotoh};
use simd::SimdBackend;

#[test]
fn scalar_and_simd_paths_match() {
    let a = b"MKTAYIAKQRQISFVKSHFSRQDILDLI";
    let b = b"MKTAYIAKQRRISFVKSHFSRQDILDLV";
    let scoring = ScoringParams::protein_default();

    let base: AlignmentResult = smith_waterman_gotoh(a, b, &scoring, SimdBackend::Scalar);
    for backend in [SimdBackend::Avx2, SimdBackend::Neon] {
        let got = smith_waterman_gotoh(a, b, &scoring, backend);
        assert_eq!(got, base, "backend mismatch: {}", backend.as_str());
    }
}
