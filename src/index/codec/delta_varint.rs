use crate::index::codec::PostingsCodec;

#[derive(Debug, Clone, Copy)]
pub struct DeltaVarintCodec;

impl PostingsCodec for DeltaVarintCodec {
    fn name(&self) -> &'static str {
        "delta_varint"
    }

    fn encode(&self, input: &[u32], out: &mut Vec<u8>) {
        let mut i = 0usize;
        while i < input.len() {
            let rem = (input.len() - i).min(4);
            out.push(rem as u8);
            for &v in &input[i..i + rem] {
                let mut x = v;
                while x >= 0x80 {
                    out.push((x as u8) | 0x80);
                    x >>= 7;
                }
                out.push(x as u8);
            }
            i += rem;
        }
    }

    fn decode(&self, bytes: &[u8], out: &mut Vec<u32>) {
        let mut i = 0usize;
        while i < bytes.len() {
            let rem = bytes[i] as usize;
            i += 1;
            for _ in 0..rem {
                let mut val = 0u32;
                let mut shift = 0u32;
                loop {
                    let b = bytes[i];
                    i += 1;
                    val |= ((b & 0x7f) as u32) << shift;
                    if b & 0x80 == 0 {
                        break;
                    }
                    shift += 7;
                }
                out.push(val);
            }
        }
    }
}
