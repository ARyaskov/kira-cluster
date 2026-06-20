use std::path::Path;

use memmap2::Mmap;

use crate::alignment::gpu::DbView;
use crate::alignment::scoring::ScoringParams;
use crate::cascade::CascadeLevel;
use crate::db::layout::DbMeta;
use crate::error::{AppError, ErrorKind, Result};
use crate::gpu::GpuBackend;
use crate::io::mmap::map_readonly;
use crate::profile::{ProfileReport, Profiler};
use crate::simd::SimdBackend;

pub mod filter;
pub mod greedy;
pub mod kmer;
pub mod output;
pub mod table;

#[derive(Debug, Clone)]
pub struct ClusterConfig {
    pub min_identity: f32,
    pub min_coverage: f32,
    /// MMseqs2-style coverage mode (0=bidirectional, 1=target, 2=query,
    /// 3=length ratio). See [`crate::cluster::filter::coverage_satisfied`].
    pub cov_mode: u8,
    pub kmer_size: Option<usize>,
    pub kmer_per_seq: usize,
    pub cascade_level: CascadeLevel,
    pub sensitivity: f32,
    pub scoring: ScoringParams,
    pub use_gpu: bool,
    pub gpu_backend: GpuBackend,
    pub batch_size: usize,
    pub gpu_memory_limit: usize,
    pub cpu_threads: usize,
    pub backend: SimdBackend,
    /// Use the reduced amino-acid alphabet for seed generation (protein only).
    /// Lets clustering find diverged homologs that exact k-mers miss.
    pub reduce_alphabet: bool,
    pub profiler: Option<std::sync::Arc<Profiler>>,
}

#[derive(Debug, Clone)]
pub struct ClusterSummary {
    pub n_clusters: usize,
    pub n_sequences: usize,
    pub profile: Option<ProfileReport>,
}

pub struct SeqDb {
    pub meta: DbMeta,
    names: Vec<String>,
    seq_offsets: Vec<u64>,
    seq_lens: Vec<u32>,
    seqs_map: Mmap,
}

impl SeqDb {
    pub fn load(dbdir: &Path) -> Result<Self> {
        let meta_map = map_readonly(&dbdir.join("meta.json"))?;
        let meta: DbMeta = serde_json::from_slice(&meta_map)
            .map_err(|e| AppError::new(ErrorKind::Parse, format!("parse meta.json: {e}")))?;

        let n = usize::try_from(meta.n_seqs)
            .map_err(|_| AppError::new(ErrorKind::Validation, "n_seqs does not fit usize"))?;

        let name_offsets = read_u64_table(&dbdir.join("name_offsets.u64"))?;
        let seq_offsets = read_u64_table(&dbdir.join("seq_offsets.u64"))?;
        let seq_lens = read_u32_table(&dbdir.join("seq_lens.u32"))?;
        if name_offsets.len() != n || seq_offsets.len() != n || seq_lens.len() != n {
            return Err(AppError::new(
                ErrorKind::Validation,
                "db tables have inconsistent length",
            ));
        }

        let names_map = map_readonly(&dbdir.join("names.bin"))?;
        let seqs_map = map_readonly(&dbdir.join("seqs.bin"))?;
        let names = decode_names(&names_map, &name_offsets)?;

        Ok(Self {
            meta,
            names,
            seq_offsets,
            seq_lens,
            seqs_map,
        })
    }

    pub fn n_seqs(&self) -> usize {
        self.seq_lens.len()
    }

    pub fn seq_len(&self, seq_id: u32) -> u32 {
        self.seq_lens[seq_id as usize]
    }

    pub fn seq(&self, seq_id: u32) -> &[u8] {
        let idx = seq_id as usize;
        let start = self.seq_offsets[idx] as usize;
        let end = start + (self.seq_lens[idx] as usize);
        &self.seqs_map[start..end]
    }

    pub fn name(&self, seq_id: u32) -> &str {
        &self.names[seq_id as usize]
    }

    pub fn default_kmer_size(&self) -> usize {
        if self.meta.dbtype == "nucleotide" {
            15
        } else {
            // Protein clustering uses the reduced (6-class) alphabet for seeds,
            // so a moderately long k retains specificity comparable to an exact
            // ~6-mer while tolerating conservative substitutions. (Exact 14-mers
            // almost never collide across diverged proteins, which is why the
            // previous default produced near-all-singletons.)
            10
        }
    }
}

impl DbView for SeqDb {
    fn seq(&self, seq_id: crate::alignment::SeqId) -> &[u8] {
        SeqDb::seq(self, seq_id)
    }
}

pub fn run_clustering(
    dbdir: &Path,
    result_dir: &Path,
    cfg: &ClusterConfig,
) -> Result<ClusterSummary> {
    let parse_guard = cfg.profiler.as_ref().map(|p| p.stage("parse_ingest"));
    let db = SeqDb::load(dbdir)?;
    drop(parse_guard);
    if let Some(p) = &cfg.profiler {
        p.set_n_sequences(db.n_seqs() as u64);
    }
    let k = cfg.kmer_size.unwrap_or_else(|| db.default_kmer_size());
    if k == 0 || cfg.kmer_per_seq == 0 {
        return Err(AppError::new(
            ErrorKind::Validation,
            "k-mer-size and kmer-per-seq must be >= 1",
        ));
    }

    let result = greedy::cluster(&db, cfg, k)?;
    let write_guard = cfg.profiler.as_ref().map(|p| p.stage("final_write"));
    output::write_outputs(result_dir, &db, &result)?;
    drop(write_guard);

    let profile = cfg.profiler.as_ref().map(|p| p.snapshot());

    Ok(ClusterSummary {
        n_clusters: result.representatives.len(),
        n_sequences: db.n_seqs(),
        profile,
    })
}

fn read_u64_table(path: &Path) -> Result<Vec<u64>> {
    let map = map_readonly(path)?;
    if map.len() % 8 != 0 {
        return Err(AppError::new(
            ErrorKind::Parse,
            format!("{} length is not divisible by 8", path.display()),
        ));
    }
    let mut out = Vec::with_capacity(map.len() / 8);
    for chunk in map.chunks_exact(8) {
        out.push(u64::from_le_bytes(chunk.try_into().expect("len 8")));
    }
    Ok(out)
}

fn read_u32_table(path: &Path) -> Result<Vec<u32>> {
    let map = map_readonly(path)?;
    if map.len() % 4 != 0 {
        return Err(AppError::new(
            ErrorKind::Parse,
            format!("{} length is not divisible by 4", path.display()),
        ));
    }
    let mut out = Vec::with_capacity(map.len() / 4);
    for chunk in map.chunks_exact(4) {
        out.push(u32::from_le_bytes(chunk.try_into().expect("len 4")));
    }
    Ok(out)
}

fn decode_names(names_map: &[u8], offsets: &[u64]) -> Result<Vec<String>> {
    let mut out = Vec::with_capacity(offsets.len());
    for (i, &off) in offsets.iter().enumerate() {
        let start = off as usize;
        let end = if i + 1 < offsets.len() {
            offsets[i + 1] as usize
        } else {
            names_map.len()
        };
        if start > end || end > names_map.len() {
            return Err(AppError::new(
                ErrorKind::Parse,
                "invalid name offsets in db",
            ));
        }
        let mut bytes = &names_map[start..end];
        if bytes.ends_with(b"\n") {
            bytes = &bytes[..bytes.len() - 1];
        }
        let s = std::str::from_utf8(bytes).map_err(|e| {
            AppError::new(ErrorKind::Parse, format!("invalid UTF-8 in names.bin: {e}"))
        })?;
        out.push(s.to_string());
    }
    Ok(out)
}
