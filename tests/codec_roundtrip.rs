#[path = "../src/simd/mod.rs"]
mod simd;

#[path = "../src/index/codec/mod.rs"]
pub mod codec;

mod index {
    pub use super::codec;
}

use index::codec::{CodecKind, delta_decode, delta_encode};

#[test]
fn codecs_roundtrip_delta_postings() {
    let ids = vec![3u32, 10, 11, 24, 100, 130, 131, 5000];
    let mut gaps = Vec::new();
    delta_encode(&ids, &mut gaps);

    for kind in [CodecKind::VByte, CodecKind::DeltaVarint, CodecKind::Bp128] {
        let mut encoded = Vec::new();
        kind.encode(&gaps, &mut encoded);

        let mut dec_gaps = Vec::new();
        kind.decode(&encoded, &mut dec_gaps);

        let mut dec_ids = Vec::new();
        delta_decode(&dec_gaps, &mut dec_ids);
        assert_eq!(dec_ids, ids, "codec={}", kind.name());
    }
}
