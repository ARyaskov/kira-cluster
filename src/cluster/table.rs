use rayon::prelude::*;

use crate::cluster::SeqDb;
use crate::cluster::kmer::extract_seeds;
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
        pool: Option<&rayon::ThreadPool>,
        profiler: Option<&Profiler>,
    ) -> Self {
        let n = db.n_seqs();
        let seed_guard = profiler.map(|p| p.stage("seed_generation"));
        let per_seq_entries: Vec<Vec<KmerEntry>> = match pool {
            Some(p) => p.install(|| {
                (0..n)
                    .into_par_iter()
                    .map(|idx| {
                        let id = idx as u32;
                        extract_seeds(db.seq(id), id, k, m)
                            .into_iter()
                            .map(|s| KmerEntry {
                                hash: s.hash,
                                seq_id: s.seq_id,
                                pos: s.pos,
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            }),
            None => (0..n)
                .map(|idx| {
                    let id = idx as u32;
                    extract_seeds(db.seq(id), id, k, m)
                        .into_iter()
                        .map(|s| KmerEntry {
                            hash: s.hash,
                            seq_id: s.seq_id,
                            pos: s.pos,
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>(),
        };
        drop(seed_guard);

        let build_guard = profiler.map(|p| p.stage("kmer_table_build"));
        let total_seeds = per_seq_entries.iter().map(|v| v.len() as u64).sum::<u64>();
        let mut entries: Vec<KmerEntry> = per_seq_entries.into_iter().flatten().collect();
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
