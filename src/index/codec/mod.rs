pub mod bp128;
pub mod delta_varint;
pub mod simd_dispatch;
pub mod vbyte;

pub trait PostingsCodec {
    fn name(&self) -> &'static str;
    fn encode(&self, input: &[u32], out: &mut Vec<u8>);
    fn decode(&self, bytes: &[u8], out: &mut Vec<u32>);
}

pub trait PostingsIterator {
    fn next(&mut self) -> Option<u32>;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CodecKind {
    VByte,
    DeltaVarint,
    Bp128,
}

impl CodecKind {
    pub fn name(self) -> &'static str {
        match self {
            CodecKind::VByte => "vbyte",
            CodecKind::DeltaVarint => "delta_varint",
            CodecKind::Bp128 => "bp128",
        }
    }

    pub fn from_name(s: &str) -> Option<Self> {
        match s {
            "vbyte" => Some(CodecKind::VByte),
            "delta_varint" => Some(CodecKind::DeltaVarint),
            "bp128" => Some(CodecKind::Bp128),
            _ => None,
        }
    }

    pub fn encode(self, input: &[u32], out: &mut Vec<u8>) {
        match self {
            CodecKind::VByte => vbyte::VByteCodec.encode(input, out),
            CodecKind::DeltaVarint => delta_varint::DeltaVarintCodec.encode(input, out),
            CodecKind::Bp128 => bp128::Bp128Codec.encode(input, out),
        }
    }

    pub fn decode(self, bytes: &[u8], out: &mut Vec<u32>) {
        match self {
            CodecKind::VByte => vbyte::VByteCodec.decode(bytes, out),
            CodecKind::DeltaVarint => delta_varint::DeltaVarintCodec.decode(bytes, out),
            CodecKind::Bp128 => bp128::Bp128Codec.decode(bytes, out),
        }
    }
}

pub fn delta_encode(sorted_ids: &[u32], out: &mut Vec<u32>) {
    out.clear();
    out.reserve(sorted_ids.len());
    let mut prev = 0u32;
    for (i, &v) in sorted_ids.iter().enumerate() {
        if i == 0 {
            out.push(v);
        } else {
            out.push(v - prev);
        }
        prev = v;
    }
}

pub fn delta_decode(gaps: &[u32], out: &mut Vec<u32>) {
    out.clear();
    out.reserve(gaps.len());
    let mut prev = 0u32;
    for (i, &g) in gaps.iter().enumerate() {
        let v = if i == 0 { g } else { prev + g };
        out.push(v);
        prev = v;
    }
}
