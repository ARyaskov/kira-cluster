use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BloomMeta {
    pub bits: u64,
    pub hashes: u32,
    pub seed: u64,
}

#[derive(Debug, Clone)]
pub struct BloomFilter {
    meta: BloomMeta,
    bits: Vec<u8>,
}

impl BloomFilter {
    pub fn build(keys: &[u64], bits_per_key: u64, hashes: u32, seed: u64) -> Self {
        let total_bits = (keys.len() as u64 * bits_per_key)
            .max(1024)
            .next_power_of_two();
        let mut bf = Self {
            meta: BloomMeta {
                bits: total_bits,
                hashes: hashes.max(1),
                seed,
            },
            bits: vec![0u8; (total_bits as usize + 7) / 8],
        };
        for &k in keys {
            bf.insert(k);
        }
        bf
    }

    pub fn from_parts(meta: BloomMeta, bits: Vec<u8>) -> Self {
        Self { meta, bits }
    }

    pub fn maybe_contains(&self, key: u64) -> bool {
        for i in 0..self.meta.hashes {
            let h = bloom_hash(key, self.meta.seed.wrapping_add(i as u64));
            let bit = (h % self.meta.bits) as usize;
            let byte = bit / 8;
            let mask = 1u8 << (bit % 8);
            if self.bits[byte] & mask == 0 {
                return false;
            }
        }
        true
    }

    pub fn meta(&self) -> &BloomMeta {
        &self.meta
    }

    pub fn bits(&self) -> &[u8] {
        &self.bits
    }

    fn insert(&mut self, key: u64) {
        for i in 0..self.meta.hashes {
            let h = bloom_hash(key, self.meta.seed.wrapping_add(i as u64));
            let bit = (h % self.meta.bits) as usize;
            let byte = bit / 8;
            let mask = 1u8 << (bit % 8);
            self.bits[byte] |= mask;
        }
    }
}

fn bloom_hash(mut x: u64, seed: u64) -> u64 {
    x ^= seed;
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ceb9fe1a85ec53);
    x ^ (x >> 33)
}
