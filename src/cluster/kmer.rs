#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KmerSeed {
    pub hash: u64,
    pub seq_id: u32,
    pub pos: u32,
}

const ROLL_BASE: u64 = 1_099_511_628_211;
const ROLL_SEED: u64 = 0x9e3779b97f4a7c15;

/// Reusable working buffers for seed extraction. Held per worker thread so the
/// per-sequence buffers are not re-allocated for every sequence (seed generation
/// dominates the fast clustering path).
#[derive(Default)]
pub struct SeedScratch {
    vals: Vec<u64>,
    invalid_prefix: Vec<u32>,
    candidates: Vec<(u64, u32)>,
}

/// Exact-alphabet seeds: every residue is its own symbol. Used by the index
/// build/search path, which relies on exact k-mer identity.
pub fn extract_seeds(seq: &[u8], seq_id: u32, k: usize, m: usize) -> Vec<KmerSeed> {
    let mut scratch = SeedScratch::default();
    let mut out = Vec::new();
    extract_into(seq, seq_id, k, m, None, &mut scratch, &mut out);
    out
}

/// Seed extraction into caller-provided, reusable buffers. `out` is cleared and
/// filled with the bottom-`m` minimizer seeds. With `reduce`, residues are
/// mapped to a 6-class Dayhoff grouping (protein clustering) so that
/// conservative substitutions collide and clustering can surface homologs well
/// below ~90% identity; otherwise every residue is its own symbol.
pub fn extract_seeds_into(
    seq: &[u8],
    seq_id: u32,
    k: usize,
    m: usize,
    reduce: bool,
    scratch: &mut SeedScratch,
    out: &mut Vec<KmerSeed>,
) {
    let table = if reduce { Some(&DAYHOFF6) } else { None };
    extract_into(seq, seq_id, k, m, table, scratch, out);
}

fn extract_into(
    seq: &[u8],
    seq_id: u32,
    k: usize,
    m: usize,
    reduce: Option<&[u8; 128]>,
    scratch: &mut SeedScratch,
    out: &mut Vec<KmerSeed>,
) {
    out.clear();
    if k == 0 || m == 0 || seq.len() < k {
        return;
    }

    let vals = &mut scratch.vals;
    let invalid_prefix = &mut scratch.invalid_prefix;
    let candidates = &mut scratch.candidates;
    vals.clear();
    invalid_prefix.clear();
    candidates.clear();

    invalid_prefix.push(0u32);
    for &b in seq {
        let up = b.to_ascii_uppercase();
        let valid = !is_ambiguous(up);
        vals.push(if valid { symbol_value(up, reduce) } else { 0 });
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

    // Keep the m smallest-hash seeds. select_nth is O(n) vs O(n log n) for a
    // full sort, and the kept SET is identical (total order on (hash, pos)); the
    // table build sorts globally afterwards, so per-sequence order is irrelevant.
    if candidates.len() > m {
        candidates.select_nth_unstable_by(m, |a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        candidates.truncate(m);
    }

    out.reserve(candidates.len());
    for &(hash, pos) in candidates.iter() {
        out.push(KmerSeed { hash, seq_id, pos });
    }
}

#[inline]
fn symbol_value(b: u8, reduce: Option<&[u8; 128]>) -> u64 {
    match reduce {
        Some(table) => (table[(b & 0x7f) as usize] as u64) + 1,
        None => (b as u64) + 1,
    }
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

/// Dayhoff 6-class reduction. Classes use distinct nonzero codes; unmapped
/// bytes share an "other" code (they are also excluded by `is_ambiguous` for
/// X/N, but B/Z/U/O are mapped to their closest class).
const DAYHOFF6: [u8; 128] = build_dayhoff6();

const fn build_dayhoff6() -> [u8; 128] {
    let mut t = [6u8; 128]; // class 6 = "other"
    // class 1: C, U (sulfur)
    t[b'C' as usize] = 1;
    t[b'U' as usize] = 1;
    // class 2: A G P S T (small)
    t[b'A' as usize] = 2;
    t[b'G' as usize] = 2;
    t[b'P' as usize] = 2;
    t[b'S' as usize] = 2;
    t[b'T' as usize] = 2;
    // class 3: D E N Q B Z (acid / amide)
    t[b'D' as usize] = 3;
    t[b'E' as usize] = 3;
    t[b'N' as usize] = 3;
    t[b'Q' as usize] = 3;
    t[b'B' as usize] = 3;
    t[b'Z' as usize] = 3;
    // class 4: H K R O (basic)
    t[b'H' as usize] = 4;
    t[b'K' as usize] = 4;
    t[b'R' as usize] = 4;
    t[b'O' as usize] = 4;
    // class 5: I L M V (hydrophobic)
    t[b'I' as usize] = 5;
    t[b'L' as usize] = 5;
    t[b'M' as usize] = 5;
    t[b'V' as usize] = 5;
    // class 7: F W Y (aromatic) — note Y/W/F kept distinct from "other"
    t[b'F' as usize] = 7;
    t[b'W' as usize] = 7;
    t[b'Y' as usize] = 7;
    t
}

#[cfg(test)]
mod tests {
    use super::{SeedScratch, extract_seeds, extract_seeds_into};

    fn reduced(seq: &[u8], id: u32) -> Vec<super::KmerSeed> {
        let mut s = SeedScratch::default();
        let mut out = Vec::new();
        extract_seeds_into(seq, id, 6, 32, true, &mut s, &mut out);
        out
    }

    #[test]
    fn reduced_alphabet_collides_conservative_substitutions() {
        // Two sequences differing only by conservative substitutions
        // (I<->L<->V, K<->R) share reduced seeds but not exact seeds.
        let a = b"MKTAYIAKQRQISFVKSHFSRQDILDLVAA";
        let b = b"MRTAYLAKQRQLSFVKSHFSRQDILDLIAA";

        let exact_a = extract_seeds(a, 0, 6, 32);
        let exact_b = extract_seeds(b, 1, 6, 32);
        let reduced_a = reduced(a, 0);
        let reduced_b = reduced(b, 1);

        let shared = |xs: &[super::KmerSeed], ys: &[super::KmerSeed]| {
            let set: std::collections::BTreeSet<u64> = xs.iter().map(|s| s.hash).collect();
            ys.iter().filter(|s| set.contains(&s.hash)).count()
        };

        let reduced_shared = shared(&reduced_a, &reduced_b);
        let exact_shared = shared(&exact_a, &exact_b);
        assert!(
            reduced_shared > exact_shared,
            "reduced should surface more shared seeds: reduced={reduced_shared} exact={exact_shared}"
        );
    }
}
