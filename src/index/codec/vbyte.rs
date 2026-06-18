use crate::index::codec::PostingsCodec;

#[derive(Debug, Clone, Copy)]
pub struct VByteCodec;

impl PostingsCodec for VByteCodec {
    fn name(&self) -> &'static str {
        "vbyte"
    }

    fn encode(&self, input: &[u32], out: &mut Vec<u8>) {
        for &mut_v in input {
            let mut v = mut_v;
            while v >= 0x80 {
                out.push((v as u8) | 0x80);
                v >>= 7;
            }
            out.push(v as u8);
        }
    }

    fn decode(&self, bytes: &[u8], out: &mut Vec<u32>) {
        let mut val = 0u32;
        let mut shift = 0u32;
        for &b in bytes {
            val |= ((b & 0x7f) as u32) << shift;
            if b & 0x80 == 0 {
                out.push(val);
                val = 0;
                shift = 0;
            } else {
                shift += 7;
            }
        }
    }
}
