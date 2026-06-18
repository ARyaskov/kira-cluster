use std::path::PathBuf;

use clap::{Args, Parser, Subcommand, ValueEnum};

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum DbTypeArg {
    Auto,
    Protein,
    Nucleotide,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum AlignmentModeArg {
    Fast,
    Sensitive,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CascadeLevelArg {
    Fast,
    Medium,
    Sensitive,
    Ultra,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum GpuBackendArg {
    Cuda,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum HotModeArg {
    Skip,
    Roaring,
    Bp128,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum SeedPolicyArg {
    FixedK,
    VariableK,
    Syncmer,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CompareMmseqsModeArg {
    Linclust,
    Cluster,
}

#[derive(Debug, Parser)]
#[command(
    name = "kira-cluster",
    version,
    about = "Deterministic clustering/search toolkit"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    #[command(name = "createdb")]
    Createdb(CreatedbArgs),
    #[command(name = "easy-cluster")]
    EasyCluster(EasyClusterArgs),
    #[command(name = "build-index")]
    BuildIndex(BuildIndexArgs),
    #[command(name = "search-index")]
    SearchIndex(SearchIndexArgs),
    #[command(name = "update-index")]
    UpdateIndex(UpdateIndexArgs),
    #[command(name = "tune-index")]
    TuneIndex(TuneIndexArgs),
    #[command(name = "index-info")]
    IndexInfo(IndexInfoArgs),
    #[command(name = "serve")]
    Serve(ServeArgs),
}

#[derive(Debug, Args, Clone)]
pub struct CommonOpts {
    #[arg(long, default_value_t = 1)]
    pub threads: usize,
    #[arg(long = "tmp-dir")]
    pub tmp_dir_opt: Option<PathBuf>,
    #[arg(long)]
    pub min_seq_id: Option<f32>,
    #[arg(short = 'c', long)]
    pub cov_mode: Option<u8>,
    #[arg(long)]
    pub cov: Option<f32>,
    #[arg(long, alias = "k-mer-per-seq")]
    pub kmer_per_seq: Option<u32>,
    #[arg(long, alias = "k-mer-size")]
    pub kmer_size: Option<u32>,
    #[arg(long)]
    pub max_seq_len: Option<u32>,
    #[arg(long, value_enum, default_value_t = DbTypeArg::Auto)]
    pub dbtype: DbTypeArg,
    #[arg(long, default_value_t = 42)]
    pub seed: u64,
    #[arg(long, default_value_t = false)]
    pub strict: bool,
    #[arg(long, default_value_t = false)]
    pub quiet: bool,
    #[arg(long, default_value_t = false)]
    pub verbose: bool,
    #[arg(long, default_value_t = false)]
    pub print_plan: bool,
    #[arg(long, default_value_t = false)]
    pub profile: bool,
    #[arg(long)]
    pub profile_json: Option<PathBuf>,
    #[arg(long, value_enum)]
    pub compare_mmseqs_mode: Option<CompareMmseqsModeArg>,

    #[arg(long, default_value_t = false)]
    pub gpu: bool,
    #[arg(long, default_value_t = false)]
    pub cuda: bool,
    #[arg(long, value_enum)]
    pub gpu_backend: Option<GpuBackendArg>,
    #[arg(long, value_enum)]
    pub alignment_mode: Option<AlignmentModeArg>,
    #[arg(long, value_enum)]
    pub cascade_level: Option<CascadeLevelArg>,
    #[arg(long)]
    pub sensitivity: Option<f32>,
    #[arg(long)]
    pub batch_size: Option<usize>,
    #[arg(long)]
    pub gpu_memory_limit: Option<usize>,
    #[arg(long)]
    pub cpu_threads: Option<usize>,
    #[arg(long)]
    pub gap_open: Option<i16>,
    #[arg(long)]
    pub gap_extend: Option<i16>,
    #[arg(long)]
    pub sub_matrix: Option<String>,
}

#[derive(Debug, Args)]
pub struct CreatedbArgs {
    pub input: PathBuf,
    pub dbdir: PathBuf,
    #[command(flatten)]
    pub opts: CommonOpts,
}

#[derive(Debug, Args)]
pub struct EasyClusterArgs {
    pub input: PathBuf,
    pub result_dir: PathBuf,
    pub tmp_dir: PathBuf,
    #[command(flatten)]
    pub opts: CommonOpts,
}

#[derive(Debug, Args)]
pub struct BuildIndexArgs {
    pub dbdir: PathBuf,
    pub indexdir: PathBuf,
    #[arg(long, default_value_t = 14)]
    pub k: usize,
    #[arg(long, default_value_t = 20)]
    pub m: usize,
    #[arg(long, default_value_t = 42)]
    pub seed: u64,
    #[arg(long)]
    pub hot_df_threshold: Option<u32>,
    #[arg(long, value_enum, default_value_t = HotModeArg::Bp128)]
    pub hot_mode: HotModeArg,
    #[arg(long, default_value_t = false)]
    pub legacy_uncompressed: bool,
    #[arg(long, default_value_t = 64)]
    pub pgm_epsilon: usize,
    #[arg(long, value_enum)]
    pub seed_policy: Option<SeedPolicyArg>,
}

#[derive(Debug, Args)]
pub struct SearchIndexArgs {
    pub indexdir: PathBuf,
    pub query_fasta: PathBuf,
    pub out_tsv: PathBuf,
    #[arg(long, default_value_t = 14)]
    pub k: usize,
    #[arg(long, default_value_t = 20)]
    pub m: usize,
    #[arg(long, default_value_t = 10)]
    pub top_k: usize,
    #[arg(long, default_value_t = 1)]
    pub min_seed_hits: u32,
    #[arg(long, default_value_t = 0.9)]
    pub min_seq_id: f32,
    #[arg(long, default_value_t = 0.8)]
    pub cov: f32,
    #[arg(long)]
    pub max_df: Option<u32>,
    #[arg(long, default_value_t = 0.90)]
    pub prune_df_quantile: f32,
    #[arg(long, default_value_t = 128)]
    pub max_seeds_per_query: u32,
    #[arg(long)]
    pub work_budget: Option<u64>,
    #[arg(long, value_enum)]
    pub seed_policy: Option<SeedPolicyArg>,
    #[arg(long)]
    pub k_mer_size_set: Option<String>,
    #[arg(long, default_value_t = 42)]
    pub seed: u64,
}

#[derive(Debug, Args)]
pub struct UpdateIndexArgs {
    pub indexdir: PathBuf,
    pub new_dbdir: PathBuf,
    #[arg(long, default_value_t = 14)]
    pub k: usize,
    #[arg(long, default_value_t = 20)]
    pub m: usize,
    #[arg(long, default_value_t = 42)]
    pub seed: u64,
    #[arg(long)]
    pub hot_df_threshold: Option<u32>,
    #[arg(long, value_enum, default_value_t = HotModeArg::Bp128)]
    pub hot_mode: HotModeArg,
    #[arg(long, default_value_t = false)]
    pub legacy_uncompressed: bool,
    #[arg(long, default_value_t = 64)]
    pub pgm_epsilon: usize,
    #[arg(long, value_enum)]
    pub seed_policy: Option<SeedPolicyArg>,
}

#[derive(Debug, Args)]
pub struct TuneIndexArgs {
    pub indexdir: PathBuf,
    #[arg(long, default_value_t = 64)]
    pub pgm_epsilon: usize,
    #[arg(long, value_enum)]
    pub seed_policy: Option<SeedPolicyArg>,
}

#[derive(Debug, Args)]
pub struct IndexInfoArgs {
    pub indexdir: PathBuf,
}

#[derive(Debug, Args)]
pub struct ServeArgs {
    pub indexdir: PathBuf,
    #[arg(long, default_value = "127.0.0.1")]
    pub host: String,
    #[arg(long, default_value_t = 8080)]
    pub port: u16,
}
