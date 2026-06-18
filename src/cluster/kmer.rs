#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KmerSeed {
    pub hash: u64,
    pub seq_id: u32,
    pub pos: u32,
}

const ROLL_BASE: u64 = 1_099_511_628_211;
const ROLL_SEED: u64 = 0x9e3779b97f4a7c15;

pub fn extract_seeds(seq: &[u8], seq_id: u32, k: usize, m: usize) -> Vec<KmerSeed> {
    if k == 0 || m == 0 || seq.len() < k {
        return Vec::new();
    }

    let mut vals = Vec::with_capacity(seq.len());
    let mut invalid_prefix = Vec::with_capacity(seq.len() + 1);
    invalid_prefix.push(0u32);
    for &b in seq {
        let up = b.to_ascii_uppercase();
        let valid = !is_ambiguous(up);
        vals.push(if valid { value(up) } else { 0 });
        let next = invalid_prefix.last().copied().unwrap_or(0) + u32::from(!valid);
        invalid_prefix.push(next);
    }

    let mut pow = 1u64;
    for _ in 1..k {
        pow = pow.wrapping_mul(ROLL_BASE);
    }

    let mut hash = ROLL_SEED;
    for &v in &vals[..k] {
        hash = hash.wrapping_mul(ROLL_BASE).wrapping_add(v);
    }

    let mut candidates: Vec<(u64, u32)> = Vec::new();
    for start in 0..=(vals.len() - k) {
        if invalid_prefix[start + k] == invalid_prefix[start] {
            candidates.push((mix_hash(hash), start as u32));
        }

        if start + k < vals.len() {
            let old = vals[start];
            let new = vals[start + k];
            hash = hash
                .wrapping_sub(old.wrapping_mul(pow))
                .wrapping_mul(ROLL_BASE)
                .wrapping_add(new);
        }
    }

    candidates.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    if candidates.len() > m {
        candidates.truncate(m);
    }

    let mut out = Vec::with_capacity(candidates.len());
    for (hash, pos) in candidates {
        out.push(KmerSeed { hash, seq_id, pos });
    }
    out
}

fn value(b: u8) -> u64 {
    (b as u64) + 1
}

fn mix_hash(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ceb9fe1a85ec53);
    x ^ (x >> 33)
}

fn is_ambiguous(b: u8) -> bool {
    !b.is_ascii_uppercase() || b == b'N' || b == b'X'
}
