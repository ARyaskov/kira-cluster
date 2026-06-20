use std::io::{Read, Write};
use std::net::TcpListener;
use std::sync::Arc;

use clap::Parser;

use crate::alignment::scoring::ScoringParams;
use crate::cascade::CascadeLevel;
use crate::cli::mmseqs_compat::{
    AlignmentModeArg, CascadeLevelArg, Cli, Commands, CommonOpts, CompareMmseqsModeArg, DbTypeArg,
    GpuBackendArg, HotModeArg, SeedPolicyArg,
};
use crate::cluster::{ClusterConfig, run_clustering};
use crate::db::createdb::{CreateDbConfig, load_meta, run_createdb};
use crate::error::{AppError, ErrorKind, Result};
use crate::gpu::{GpuBackend, is_backend_available};
use crate::index::IndexHandle;
use crate::index::build::{BuildIndexConfig, build_fresh_index};
use crate::index::index_info::print_index_info;
use crate::index::search::{SearchOpts, search_index};
use crate::index::tier::HotMode;
use crate::index::tune::seed_policy::SeedPolicyMode;
use crate::index::tune::{TuneConfig, tune_index};
use crate::index::update::update_index;
use crate::profile::Profiler;
use crate::seq::DbType;
use crate::simd;
use crate::util::log::{self, LogConfig};
use crate::util::paths::{is_db_dir, query_db_dir};

pub mod mmseqs_compat;

const WARN_UNIMPLEMENTED: &str = "KW2001";

fn dbtype_from_arg(v: DbTypeArg) -> DbType {
    match v {
        DbTypeArg::Auto => DbType::Auto,
        DbTypeArg::Protein => DbType::Protein,
        DbTypeArg::Nucleotide => DbType::Nucleotide,
    }
}

fn validate_common(opts: &CommonOpts) -> Result<Vec<String>> {
    if opts.threads == 0 {
        return Err(AppError::new(
            ErrorKind::Validation,
            "--threads must be >= 1",
        ));
    }
    if let Some(v) = opts.min_seq_id {
        if !(0.0..=1.0).contains(&v) {
            return Err(AppError::new(
                ErrorKind::Validation,
                "--min-seq-id must be within 0..1",
            ));
        }
    }
    if let Some(v) = opts.cov_mode {
        if v > 3 {
            return Err(AppError::new(
                ErrorKind::Validation,
                "--cov-mode must be one of 0,1,2,3",
            ));
        }
    }
    if let Some(v) = opts.cov {
        if !(0.0..=1.0).contains(&v) {
            return Err(AppError::new(
                ErrorKind::Validation,
                "--cov must be within 0..1",
            ));
        }
    }
    if let Some(v) = opts.kmer_size {
        if v == 0 {
            return Err(AppError::new(
                ErrorKind::Validation,
                "--kmer-size must be >= 1",
            ));
        }
    }
    if let Some(v) = opts.kmer_per_seq {
        if v == 0 {
            return Err(AppError::new(
                ErrorKind::Validation,
                "--kmer-per-seq must be >= 1",
            ));
        }
    }
    if let Some(v) = opts.sensitivity {
        if v < 0.0 {
            return Err(AppError::new(
                ErrorKind::Validation,
                "--sensitivity must be >= 0",
            ));
        }
    }
    if let Some(v) = opts.batch_size {
        if v == 0 {
            return Err(AppError::new(
                ErrorKind::Validation,
                "--batch-size must be >= 1",
            ));
        }
    }
    if let Some(v) = opts.cpu_threads {
        if v == 0 {
            return Err(AppError::new(
                ErrorKind::Validation,
                "--cpu-threads must be >= 1",
            ));
        }
    }

    let mut warnings = Vec::new();
    if opts.gpu || opts.cuda {
        warnings.push("GPU_FLAG_REQUESTED".to_string());
    }
    if let Some(m) = &opts.sub_matrix {
        let ml = m.to_ascii_lowercase();
        if ml != "simple" && ml != "blosum62" {
            warnings.push(format!("UNIMPLEMENTED_SUB_MATRIX {m}"));
        }
    }
    warnings.sort();
    Ok(warnings)
}

fn emit_warnings(cfg: LogConfig, warnings: &[String], strict: bool) -> Result<()> {
    for flag in warnings {
        if flag.starts_with("UNIMPLEMENTED_") || flag.starts_with("GPU_") {
            log::warn(cfg, WARN_UNIMPLEMENTED, flag);
        } else {
            log::warn(
                cfg,
                WARN_UNIMPLEMENTED,
                &format!("UNIMPLEMENTED_FLAG {flag}"),
            );
        }
    }
    if strict && !warnings.is_empty() {
        return Err(AppError::new(
            ErrorKind::Unsupported,
            "strict mode: unimplemented flags present",
        ));
    }
    Ok(())
}

fn print_common_plan(opts: &CommonOpts) {
    println!("plan.threads={}", opts.threads);
    println!(
        "plan.tmp_dir_opt={}",
        opts.tmp_dir_opt
            .as_ref()
            .map_or("none".to_string(), |p| p.display().to_string())
    );
    println!(
        "plan.min_seq_id={}",
        opts.min_seq_id
            .map_or("none".to_string(), |v| v.to_string())
    );
    println!(
        "plan.cov_mode={}",
        opts.cov_mode.map_or("none".to_string(), |v| v.to_string())
    );
    println!(
        "plan.cov={}",
        opts.cov.map_or("none".to_string(), |v| v.to_string())
    );
    println!(
        "plan.kmer_per_seq={}",
        opts.kmer_per_seq
            .map_or("none".to_string(), |v| v.to_string())
    );
    println!(
        "plan.kmer_size={}",
        opts.kmer_size.map_or("none".to_string(), |v| v.to_string())
    );
    println!(
        "plan.max_seq_len={}",
        opts.max_seq_len
            .map_or("none".to_string(), |v| v.to_string())
    );
    println!(
        "plan.requested_dbtype={}",
        dbtype_from_arg(opts.dbtype).as_str()
    );
    println!("plan.seed={}", opts.seed);
    println!("plan.strict={}", opts.strict);
    println!("plan.quiet={}", opts.quiet);
    println!("plan.verbose={}", opts.verbose);
    println!("plan.profile={}", opts.profile);
    println!(
        "plan.profile_json={}",
        opts.profile_json
            .as_ref()
            .map_or("none".to_string(), |p| p.display().to_string())
    );
    println!(
        "plan.compare_mmseqs_mode={}",
        opts.compare_mmseqs_mode
            .map_or("none".to_string(), |m| match m {
                CompareMmseqsModeArg::Linclust => "linclust".to_string(),
                CompareMmseqsModeArg::Cluster => "cluster".to_string(),
            })
    );
    println!("plan.gpu={}", opts.gpu);
    println!("plan.cuda={}", opts.cuda);
    println!(
        "plan.gpu_backend={}",
        opts.gpu_backend.map_or("none".to_string(), |v| match v {
            GpuBackendArg::Cuda => "cuda".to_string(),
        })
    );
    println!(
        "plan.alignment_mode={}",
        opts.alignment_mode.map_or("none".to_string(), |v| match v {
            AlignmentModeArg::Fast => "fast".to_string(),
            AlignmentModeArg::Sensitive => "sensitive".to_string(),
        })
    );
    println!(
        "plan.sensitivity={}",
        opts.sensitivity
            .map_or("none".to_string(), |v| v.to_string())
    );
    println!(
        "plan.cascade_level={}",
        opts.cascade_level.map_or("none".to_string(), |v| match v {
            CascadeLevelArg::Fast => "fast".to_string(),
            CascadeLevelArg::Medium => "medium".to_string(),
            CascadeLevelArg::Sensitive => "sensitive".to_string(),
            CascadeLevelArg::Ultra => "ultra".to_string(),
        })
    );
    println!(
        "plan.batch_size={}",
        opts.batch_size
            .map_or("none".to_string(), |v| v.to_string())
    );
    println!(
        "plan.gpu_memory_limit={}",
        opts.gpu_memory_limit
            .map_or("none".to_string(), |v| v.to_string())
    );
    println!(
        "plan.cpu_threads={}",
        opts.cpu_threads
            .map_or("none".to_string(), |v| v.to_string())
    );
    println!(
        "plan.gap_open={}",
        opts.gap_open.map_or("none".to_string(), |v| v.to_string())
    );
    println!(
        "plan.gap_extend={}",
        opts.gap_extend
            .map_or("none".to_string(), |v| v.to_string())
    );
    println!(
        "plan.sub_matrix={}",
        opts.sub_matrix
            .clone()
            .unwrap_or_else(|| "none".to_string())
    );
}

fn print_mmseqs_comparison_header(opts: &CommonOpts) {
    if let Some(mode) = opts.compare_mmseqs_mode {
        println!("Comparison parameters:");
        println!(
            "  mode: {}",
            match mode {
                CompareMmseqsModeArg::Linclust => "linclust",
                CompareMmseqsModeArg::Cluster => "cluster",
            }
        );
        println!("  min_seq_id: {}", opts.min_seq_id.unwrap_or(0.9));
        println!("  cov_mode: {}", opts.cov_mode.unwrap_or(0));
        println!("  coverage: {}", opts.cov.unwrap_or(0.8));
    }
}

fn resolve_cascade_level(opts: &CommonOpts) -> CascadeLevel {
    if let Some(v) = opts.cascade_level {
        return match v {
            CascadeLevelArg::Fast => CascadeLevel::Fast,
            CascadeLevelArg::Medium => CascadeLevel::Medium,
            CascadeLevelArg::Sensitive => CascadeLevel::Sensitive,
            CascadeLevelArg::Ultra => CascadeLevel::Ultra,
        };
    }
    // Default to the fast (positional) cascade for predictable speed on any
    // input size. The reduced-alphabet candidate generation and BLOSUM62 scoring
    // still apply; pass --alignment-mode sensitive / --cascade-level sensitive to
    // enable gapped alignment (catches indel-containing homologs, slower).
    match opts.alignment_mode {
        Some(AlignmentModeArg::Fast) => CascadeLevel::Fast,
        Some(AlignmentModeArg::Sensitive) => CascadeLevel::Sensitive,
        None => CascadeLevel::Fast,
    }
}

fn resolve_gpu_backend(opts: &CommonOpts) -> GpuBackend {
    match opts.gpu_backend.unwrap_or(GpuBackendArg::Cuda) {
        GpuBackendArg::Cuda => GpuBackend::Cuda,
    }
}

fn resolve_scoring(opts: &CommonOpts, resolved_dbtype: &str) -> ScoringParams {
    let mut scoring = if resolved_dbtype == "nucleotide" {
        ScoringParams::nucleotide_default()
    } else {
        // Protein defaults to BLOSUM62; `--sub-matrix simple` selects the flat
        // match/mismatch model.
        match opts.sub_matrix.as_deref().map(str::to_ascii_lowercase) {
            Some(ref s) if s == "simple" => ScoringParams::protein_default(),
            _ => ScoringParams::protein_blosum62(),
        }
    };
    if let Some(v) = opts.gap_open {
        scoring.gap_open = v;
    }
    if let Some(v) = opts.gap_extend {
        scoring.gap_extend = v;
    }
    scoring
}

fn resolve_hot_mode(v: HotModeArg) -> HotMode {
    match v {
        HotModeArg::Skip => HotMode::Skip,
        HotModeArg::Roaring => HotMode::Roaring,
        HotModeArg::Bp128 => HotMode::Bp128,
    }
}

fn resolve_seed_policy_mode(v: Option<SeedPolicyArg>) -> Option<SeedPolicyMode> {
    v.map(|m| match m {
        SeedPolicyArg::FixedK => SeedPolicyMode::FixedK,
        SeedPolicyArg::VariableK => SeedPolicyMode::VariableK,
        SeedPolicyArg::Syncmer => SeedPolicyMode::Syncmer,
    })
}

pub fn run() -> Result<()> {
    let cli = Cli::parse();
    let backend = simd::active_backend();

    match cli.command {
        Commands::Createdb(args) => {
            let log_cfg = LogConfig::new(args.opts.quiet, args.opts.verbose);
            log::backend_startup(log_cfg, backend);

            let warnings = validate_common(&args.opts)?;
            emit_warnings(log_cfg, &warnings, args.opts.strict)?;

            let meta = run_createdb(&CreateDbConfig {
                input: args.input,
                dbdir: args.dbdir,
                dbtype: dbtype_from_arg(args.opts.dbtype),
                max_seq_len: args.opts.max_seq_len,
            })?;

            if args.opts.print_plan {
                println!("plan.command=createdb");
                println!("plan.backend={}", backend.as_str());
                println!("plan.resolved_dbtype={}", meta.dbtype);
                println!("plan.n_seqs={}", meta.n_seqs);
                println!("plan.total_bases={}", meta.total_bases);
                print_common_plan(&args.opts);
            }

            Ok(())
        }
        Commands::EasyCluster(args) => {
            let log_cfg = LogConfig::new(args.opts.quiet, args.opts.verbose);
            log::backend_startup(log_cfg, backend);

            let warnings = validate_common(&args.opts)?;
            let mut effective_warnings = warnings;
            let use_gpu = args.opts.gpu || args.opts.cuda;
            let gpu_backend = resolve_gpu_backend(&args.opts);
            if use_gpu && !is_backend_available(gpu_backend) {
                effective_warnings.push(format!(
                    "GPU_BACKEND_UNAVAILABLE_FALLBACK {}",
                    gpu_backend.as_str()
                ));
            }
            emit_warnings(log_cfg, &effective_warnings, args.opts.strict)?;

            std::fs::create_dir_all(&args.result_dir).map_err(|e| {
                AppError::io(format!("create dir {}", args.result_dir.display()), e)
            })?;

            let effective_tmp = args
                .opts
                .tmp_dir_opt
                .clone()
                .unwrap_or_else(|| args.tmp_dir.clone());
            std::fs::create_dir_all(&effective_tmp)
                .map_err(|e| AppError::io(format!("create dir {}", effective_tmp.display()), e))?;

            let query_db = if is_db_dir(&args.input) {
                args.input.clone()
            } else {
                let dbdir = query_db_dir(&effective_tmp);
                let _meta = run_createdb(&CreateDbConfig {
                    input: args.input.clone(),
                    dbdir: dbdir.clone(),
                    dbtype: dbtype_from_arg(args.opts.dbtype),
                    max_seq_len: args.opts.max_seq_len,
                })?;
                dbdir
            };

            let meta = load_meta(&query_db)?;
            let profiler = if args.opts.profile || args.opts.profile_json.is_some() {
                Some(Arc::new(Profiler::new(true)))
            } else {
                None
            };

            if args.opts.print_plan {
                println!("plan.command=easy-cluster");
                println!("plan.backend={}", backend.as_str());
                println!("plan.input_db={}", query_db.display());
                println!("plan.result_dir={}", args.result_dir.display());
                println!("plan.tmp_dir={}", effective_tmp.display());
                println!("plan.resolved_dbtype={}", meta.dbtype);
                println!("plan.n_seqs={}", meta.n_seqs);
                print_common_plan(&args.opts);
                println!("plan.stage=cluster_v1_enabled");
            }
            print_mmseqs_comparison_header(&args.opts);

            let summary = run_clustering(
                &query_db,
                &args.result_dir,
                &ClusterConfig {
                    min_identity: args.opts.min_seq_id.unwrap_or(0.9),
                    min_coverage: args.opts.cov.unwrap_or(0.8),
                    cov_mode: args.opts.cov_mode.unwrap_or(0),
                    kmer_size: args.opts.kmer_size.map(|v| v as usize),
                    kmer_per_seq: ((args.opts.kmer_per_seq.unwrap_or(20) as f32)
                        * args.opts.sensitivity.unwrap_or(1.0).max(1.0))
                    .round() as usize,
                    cascade_level: resolve_cascade_level(&args.opts),
                    sensitivity: args.opts.sensitivity.unwrap_or(1.0),
                    scoring: resolve_scoring(&args.opts, &meta.dbtype),
                    use_gpu: use_gpu && is_backend_available(gpu_backend),
                    gpu_backend,
                    batch_size: args.opts.batch_size.unwrap_or(1024),
                    gpu_memory_limit: args.opts.gpu_memory_limit.unwrap_or(1 << 30),
                    cpu_threads: args.opts.cpu_threads.unwrap_or(args.opts.threads),
                    backend,
                    reduce_alphabet: meta.dbtype == "protein",
                    profiler: profiler.clone(),
                },
            )?;

            if !args.opts.quiet {
                eprintln!(
                    "INFO CLUSTER_DONE clusters={} sequences={}",
                    summary.n_clusters, summary.n_sequences
                );
            }

            if let Some(report) = &summary.profile {
                if args.opts.profile {
                    report.print_console();
                }
                if let Some(path) = &args.opts.profile_json {
                    report.write_json(path)?;
                }
            }

            Ok(())
        }
        Commands::BuildIndex(args) => {
            build_fresh_index(&BuildIndexConfig {
                dbdir: args.dbdir,
                indexdir: args.indexdir,
                k: args.k,
                m: args.m,
                seed: args.seed,
                hot_df_threshold: args.hot_df_threshold,
                hot_mode: resolve_hot_mode(args.hot_mode),
                legacy_uncompressed: args.legacy_uncompressed,
                pgm_epsilon: args.pgm_epsilon,
                seed_policy_mode: resolve_seed_policy_mode(args.seed_policy),
            })?;
            Ok(())
        }
        Commands::SearchIndex(args) => {
            let handle = IndexHandle::open(&args.indexdir)?;
            if matches!(args.seed_policy, Some(SeedPolicyArg::Syncmer)) {
                eprintln!(
                    "WARN {} UNIMPLEMENTED_FLAG --seed-policy syncmer",
                    WARN_UNIMPLEMENTED
                );
            }
            let scoring = match handle
                .segments
                .first()
                .map(|s| s.db.meta.dbtype.as_str())
            {
                Some("nucleotide") => ScoringParams::nucleotide_default(),
                _ => ScoringParams::protein_blosum62(),
            };
            let opts = SearchOpts {
                k: args.k,
                m: args.m,
                top_k: args.top_k,
                min_seed_hits: args.min_seed_hits,
                min_identity: args.min_seq_id,
                min_coverage: args.cov,
                max_df: args.max_df,
                prune_df_quantile: args.prune_df_quantile,
                max_seeds_per_query: args.max_seeds_per_query,
                work_budget: args.work_budget,
                scoring,
                backend,
                seed: args.seed,
                seed_policy: resolve_seed_policy_mode(args.seed_policy)
                    .map(|mode| crate::index::tune::seed_policy::SeedPolicy {
                        mode,
                        k_values: vec![
                            args.k.saturating_sub(1) as u32,
                            args.k as u32,
                            (args.k + 1) as u32,
                        ],
                        m: args.m as u32,
                        rules: crate::index::tune::seed_policy::SeedPolicyRules {
                            short_query_use_k: args.k.saturating_sub(1) as u32,
                            long_query_use_k: args.k as u32,
                        },
                    })
                    .or_else(|| handle.seed_policy.clone()),
            };
            let stats = search_index(&handle, &args.query_fasta, &args.out_tsv, &opts)?;
            eprintln!(
                "INFO SEARCH_DONE queries={} hits={} zero_hit_queries={}",
                stats.queries,
                stats.total_hits,
                stats.zero_hit_queries()
            );
            Ok(())
        }
        Commands::UpdateIndex(args) => {
            update_index(
                &args.indexdir,
                BuildIndexConfig {
                    dbdir: args.new_dbdir,
                    indexdir: args.indexdir.clone(),
                    k: args.k,
                    m: args.m,
                    seed: args.seed,
                    hot_df_threshold: args.hot_df_threshold,
                    hot_mode: resolve_hot_mode(args.hot_mode),
                    legacy_uncompressed: args.legacy_uncompressed,
                    pgm_epsilon: args.pgm_epsilon,
                    seed_policy_mode: resolve_seed_policy_mode(args.seed_policy),
                },
            )?;
            Ok(())
        }
        Commands::TuneIndex(args) => {
            if matches!(args.seed_policy, Some(SeedPolicyArg::Syncmer)) {
                eprintln!(
                    "WARN {} UNIMPLEMENTED_FLAG --seed-policy syncmer",
                    WARN_UNIMPLEMENTED
                );
            }
            tune_index(&TuneConfig {
                indexdir: args.indexdir,
                pgm_epsilon: args.pgm_epsilon,
                seed_policy_mode: resolve_seed_policy_mode(args.seed_policy).map(|m| match m {
                    SeedPolicyMode::FixedK => "fixed_k".to_string(),
                    SeedPolicyMode::VariableK => "variable_k".to_string(),
                    SeedPolicyMode::Syncmer => "syncmer".to_string(),
                }),
            })
        }
        Commands::IndexInfo(args) => {
            print_index_info(&args.indexdir)?;
            Ok(())
        }
        Commands::Serve(args) => {
            let _handle = IndexHandle::open(&args.indexdir)?;
            eprintln!(
                "WARN {} EXPERIMENTAL serve is a stub: it answers a single request with a static JSON body and exits",
                WARN_UNIMPLEMENTED
            );
            let addr = format!("{}:{}", args.host, args.port);
            let listener =
                TcpListener::bind(&addr).map_err(|e| AppError::io(format!("bind {}", addr), e))?;
            for stream in listener.incoming().take(1) {
                let mut stream = stream.map_err(|e| AppError::io("accept", e))?;
                let mut buf = [0u8; 1024];
                let _ = stream.read(&mut buf);
                let body = b"{\"status\":\"ok\",\"message\":\"serve stub\"}";
                let resp = format!(
                    "HTTP/1.1 200 OK\r\nContent-Type: application/json\r\nContent-Length: {}\r\n\r\n",
                    body.len()
                );
                stream
                    .write_all(resp.as_bytes())
                    .and_then(|_| stream.write_all(body))
                    .map_err(|e| AppError::io("write response", e))?;
            }
            Ok(())
        }
    }
}
