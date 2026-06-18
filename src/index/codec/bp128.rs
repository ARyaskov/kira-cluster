use crate::index::codec::PostingsCodec;

#[derive(Debug, Clone, Copy)]
pub struct Bp128Codec;

impl PostingsCodec for Bp128Codec {
    fn name(&self) -> &'static str {
        "bp128"
    }

    fn encode(&self, input: &[u32], out: &mut Vec<u8>) {
        let mut i = 0usize;
        while i < input.len() {
            let n = (input.len() - i).min(128);
            out.push(n as u8);
            let block = &input[i..i + n];
            let max_v = block.iter().copied().max().unwrap_or(0);
            let bits = (32 - max_v.leading_zeros()).max(1) as u8;
            out.push(bits);
            pack_block(block, bits, out);
            i += n;
        }
    }

    fn decode(&self, bytes: &[u8], out: &mut Vec<u32>) {
        let mut i = 0usize;
        while i < bytes.len() {
            let n = bytes[i] as usize;
            let bits = bytes[i + 1];
            i += 2;
            let words = ((n * bits as usize) + 31) / 32;
            let words_bytes = words * 4;
            unpack_block(&bytes[i..i + words_bytes], n, bits, out);
            i += words_bytes;
        }
    }
}

fn pack_block(input: &[u32], bits: u8, out: &mut Vec<u8>) {
    let mut bitbuf: u64 = 0;
    let mut bitcount = 0u32;
    let mask = if bits == 32 {
        u64::MAX
    } else {
        (1u64 << bits) - 1
    };

    for &v in input {
        bitbuf |= ((v as u64) & mask) << bitcount;
        bitcount += bits as u32;
        while bitcount >= 32 {
            out.extend_from_slice(&(bitbuf as u32).to_le_bytes());
            bitbuf >>= 32;
            bitcount -= 32;
        }
    }

    if bitcount > 0 {
        out.extend_from_slice(&(bitbuf as u32).to_le_bytes());
    }
}

fn unpack_block(bytes: &[u8], n: usize, bits: u8, out: &mut Vec<u32>) {
    let mut words = Vec::with_capacity(bytes.len() / 4);
    for c in bytes.chunks_exact(4) {
        words.push(u32::from_le_bytes(c.try_into().expect("4")) as u64);
    }

    let mut word_idx = 0usize;
    let mut bitbuf = if words.is_empty() { 0 } else { words[0] };
    let mut bits_in_buf = 32u32;
    let mask = if bits == 32 {
        u64::MAX
    } else {
        (1u64 << bits) - 1
    };

    for _ in 0..n {
        while bits_in_buf < bits as u32 {
            word_idx += 1;
            let next = if word_idx < words.len() {
                words[word_idx]
            } else {
                0
            };
            bitbuf |= next << bits_in_buf;
            bits_in_buf += 32;
        }
        out.push((bitbuf & mask) as u32);
        bitbuf >>= bits;
        bits_in_buf -= bits as u32;
    }
}
