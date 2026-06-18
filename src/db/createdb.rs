use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::db::layout::DbMeta;
use crate::error::{AppError, ErrorKind, Result};
use crate::io::atomic::write_atomic;
use crate::io::mmap::map_readonly;
use crate::seq::{DbType, InputFormat, detect_input_format, fasta, fastq, validate_seq};
use crate::simd;

#[derive(Debug, Clone)]
pub struct CreateDbConfig {
    pub input: PathBuf,
    pub dbdir: PathBuf,
    pub dbtype: DbType,
    pub max_seq_len: Option<u32>,
}

struct TempWriters {
    names_tmp: PathBuf,
    seqs_tmp: PathBuf,
    names: BufWriter<File>,
    seqs: BufWriter<File>,
}

fn open_temp_writers(dbdir: &Path) -> Result<TempWriters> {
    let names_tmp = dbdir.join("names.bin.tmp");
    let seqs_tmp = dbdir.join("seqs.bin.tmp");
    let names = BufWriter::new(
        File::create(&names_tmp)
            .map_err(|e| AppError::io(format!("create {}", names_tmp.display()), e))?,
    );
    let seqs = BufWriter::new(
        File::create(&seqs_tmp)
            .map_err(|e| AppError::io(format!("create {}", seqs_tmp.display()), e))?,
    );

    Ok(TempWriters {
        names_tmp,
        seqs_tmp,
        names,
        seqs,
    })
}

fn finalize_temp_file(tmp: &Path, final_path: &Path) -> Result<()> {
    fs::rename(tmp, final_path).map_err(|e| {
        AppError::io(
            format!("rename {} -> {}", tmp.display(), final_path.display()),
            e,
        )
    })
}

fn u64_vec_to_le_bytes(v: &[u64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(v.len() * 8);
    for x in v {
        out.extend_from_slice(&x.to_le_bytes());
    }
    out
}

fn u32_vec_to_le_bytes(v: &[u32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(v.len() * 4);
    for x in v {
        out.extend_from_slice(&x.to_le_bytes());
    }
    out
}

pub fn run_createdb(cfg: &CreateDbConfig) -> Result<DbMeta> {
    fs::create_dir_all(&cfg.dbdir)
        .map_err(|e| AppError::io(format!("create dir {}", cfg.dbdir.display()), e))?;

    let mut writers = open_temp_writers(&cfg.dbdir)?;

    let mut n_seqs: u64 = 0;
    let mut total_bases: u64 = 0;
    let mut observed_max_len: u32 = 0;
    let mut names_pos: u64 = 0;
    let mut seqs_pos: u64 = 0;
    let mut name_offsets: Vec<u64> = Vec::new();
    let mut seq_offsets: Vec<u64> = Vec::new();
    let mut seq_lens: Vec<u32> = Vec::new();
    let mut checksums: Vec<u64> = Vec::new();
    let mut all_nucleotide = true;

    let format = detect_input_format(&cfg.input)?;

    let mut on_record = |name: String, seq: Vec<u8>| -> Result<()> {
        let len_u32 = u32::try_from(seq.len())
            .map_err(|_| AppError::new(ErrorKind::Validation, "sequence length exceeds u32"))?;

        if let Some(limit) = cfg.max_seq_len {
            if len_u32 > limit {
                return Err(AppError::new(
                    ErrorKind::Validation,
                    format!("sequence '{name}' exceeds max length {limit}"),
                ));
            }
        }

        validate_seq(&seq, cfg.dbtype)?;

        if all_nucleotide {
            const ALLOWED: &[u8] = b"ACGTUNRYSWKMBDHV";
            all_nucleotide = seq.iter().all(|b| ALLOWED.contains(b));
        }

        let checksum = simd::hash_bytes32(&seq);

        name_offsets.push(names_pos);
        seq_offsets.push(seqs_pos);
        seq_lens.push(len_u32);
        checksums.push(checksum);

        writers
            .names
            .write_all(name.as_bytes())
            .and_then(|_| writers.names.write_all(b"\n"))
            .map_err(|e| AppError::io("write names.bin.tmp", e))?;
        writers
            .seqs
            .write_all(&seq)
            .map_err(|e| AppError::io("write seqs.bin.tmp", e))?;

        names_pos += (name.len() as u64) + 1;
        seqs_pos += seq.len() as u64;
        n_seqs += 1;
        total_bases += seq.len() as u64;
        observed_max_len = observed_max_len.max(len_u32);

        Ok(())
    };

    match format {
        InputFormat::Fasta => fasta::parse_fasta(&cfg.input, &mut on_record)?,
        InputFormat::Fastq => fastq::parse_fastq(&cfg.input, &mut on_record)?,
    }

    writers
        .names
        .flush()
        .map_err(|e| AppError::io("flush names.bin.tmp", e))?;
    writers
        .seqs
        .flush()
        .map_err(|e| AppError::io("flush seqs.bin.tmp", e))?;

    finalize_temp_file(&writers.names_tmp, &cfg.dbdir.join("names.bin"))?;
    finalize_temp_file(&writers.seqs_tmp, &cfg.dbdir.join("seqs.bin"))?;

    write_atomic(
        &cfg.dbdir.join("name_offsets.u64"),
        &u64_vec_to_le_bytes(&name_offsets),
    )?;
    write_atomic(
        &cfg.dbdir.join("seq_offsets.u64"),
        &u64_vec_to_le_bytes(&seq_offsets),
    )?;
    write_atomic(
        &cfg.dbdir.join("seq_lens.u32"),
        &u32_vec_to_le_bytes(&seq_lens),
    )?;
    write_atomic(
        &cfg.dbdir.join("checksums.u64"),
        &u64_vec_to_le_bytes(&checksums),
    )?;

    let resolved_dbtype = match cfg.dbtype {
        DbType::Auto => {
            if all_nucleotide {
                "nucleotide".to_string()
            } else {
                "protein".to_string()
            }
        }
        DbType::Protein => "protein".to_string(),
        DbType::Nucleotide => "nucleotide".to_string(),
    };

    let meta = DbMeta {
        version: 1,
        dbtype: resolved_dbtype,
        n_seqs,
        total_bases,
        max_seq_len: observed_max_len,
        encoding: "ascii_upper".to_string(),
        endian: "little".to_string(),
    };

    let meta_bytes = serde_json::to_vec_pretty(&meta)
        .map_err(|e| AppError::new(ErrorKind::Internal, format!("serialize meta.json: {e}")))?;
    write_atomic(&cfg.dbdir.join("meta.json"), &meta_bytes)?;

    Ok(meta)
}

pub fn load_meta(dbdir: &Path) -> Result<DbMeta> {
    let map = map_readonly(&dbdir.join("meta.json"))?;
    let meta: DbMeta = serde_json::from_slice(&map)
        .map_err(|e| AppError::new(ErrorKind::Parse, format!("parse meta.json: {e}")))?;
    Ok(meta)
}
