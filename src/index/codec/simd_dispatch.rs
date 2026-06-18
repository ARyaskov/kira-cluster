use crate::index::codec::bp128::Bp128Codec;
use crate::index::codec::{CodecKind, PostingsCodec};
use crate::simd::SimdBackend;

pub fn decode_bp128_simd(bytes: &[u8], out: &mut Vec<u32>, _backend: SimdBackend) {
    Bp128Codec.decode(bytes, out);
}

pub fn preferred_codec_for_cold() -> CodecKind {
    CodecKind::Bp128
}
