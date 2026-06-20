use rayon::prelude::*;

use crate::cluster::SeqDb;
use crate::cluster::kmer::{KmerSeed, SeedScratch, extract_seeds_into};
use crate::profile::{Profiler, WorkCounters};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KmerEntry {
    pub hash: u64,
    pub seq_id: u32,
    pub pos: u32,
}

#[derive(Debug, Clone)]
pub struct KmerTable {
    entries: Vec<KmerEntry>,
}

impl KmerTable {
    pub fn build(
        db: &SeqDb,
        k: usize,
        m: usize,
        reduce: bool,
        pool: Option<&rayon::ThreadPool>,
        profiler: Option<&Profiler>,
    ) -> Self {
        let n = db.n_seqs();
        let to_entries = |buf: &[KmerSeed]| -> Vec<KmerEntry> {
            buf.iter()
                .map(|s| KmerEntry {
                    hash: s.hash,
                    seq_id: s.seq_id,
                    pos: s.pos,
                })
                .collect()
        };
        // Reuse per-thread scratch buffers across sequences (seed generation
        // dominates the fast path); order is preserved, so output stays
        // deterministic.
        let seed_guard = profiler.map(|p| p.stage("seed_generation"));
        let mut entries: Vec<KmerEntry> = match pool {
            Some(p) => p.install(|| {
                (0..n)
                    .into_par_iter()
                    .map_init(
                        || (SeedScratch::default(), Vec::<KmerSeed>::new()),
                        |(scratch, buf), idx| {
                            extract_seeds_into(db.seq(idx as u32), idx as u32, k, m, reduce, scratch, buf);
                            to_entries(buf)
                        },
                    )
                    .flatten()
                    .collect()
            }),
            None => {
                let mut scratch = SeedScratch::default();
                let mut buf = Vec::new();
                let mut acc = Vec::new();
                for idx in 0..n {
                    extract_seeds_into(db.seq(idx as u32), idx as u32, k, m, reduce, &mut scratch, &mut buf);
                    acc.extend(buf.iter().map(|s| KmerEntry {
                        hash: s.hash,
                        seq_id: s.seq_id,
                        pos: s.pos,
                    }));
                }
                acc
            }
        };
        drop(seed_guard);

        let build_guard = profiler.map(|p| p.stage("kmer_table_build"));
        let total_seeds = entries.len() as u64;
        drop(build_guard);

        if let Some(p) = profiler {
            p.add_counters(&WorkCounters {
                total_seeds_emitted: total_seeds,
                ..WorkCounters::default()
            });
        }

        let sort_guard = profiler.map(|p| p.stage("kmer_sort"));
        entries.sort_unstable_by(|a, b| {
            a.hash
                .cmp(&b.hash)
                .then(a.seq_id.cmp(&b.seq_id))
                .then(a.pos.cmp(&b.pos))
        });
        drop(sort_guard);

        if let Some(p) = profiler {
            let mut unique = 0u64;
            let mut prev = None;
            for e in &entries {
                if prev != Some(e.hash) {
                    unique += 1;
                    prev = Some(e.hash);
                }
            }
            p.set_unique_keys(unique);
        }

        Self { entries }
    }

    pub fn groups(&self) -> KmerGroupIter<'_> {
        KmerGroupIter {
            entries: &self.entries,
            idx: 0,
        }
    }
}

pub struct KmerGroupIter<'a> {
    entries: &'a [KmerEntry],
    idx: usize,
}

impl<'a> Iterator for KmerGroupIter<'a> {
    type Item = &'a [KmerEntry];

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.entries.len() {
            return None;
        }
        let start = self.idx;
        let hash = self.entries[start].hash;
        while self.idx < self.entries.len() && self.entries[self.idx].hash == hash {
            self.idx += 1;
        }
        Some(&self.entries[start..self.idx])
    }
}
