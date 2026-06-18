use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::error::{AppError, ErrorKind, Result};

pub mod fasta;
pub mod fastq;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputFormat {
    Fasta,
    Fastq,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DbType {
    Auto,
    Protein,
    Nucleotide,
}

impl DbType {
    pub fn as_str(self) -> &'static str {
        match self {
            DbType::Auto => "auto",
            DbType::Protein => "protein",
            DbType::Nucleotide => "nucleotide",
        }
    }
}

#[derive(Debug, Clone)]
pub struct SeqRecord {
    pub id: u32,
    pub name: String,
    pub seq: Vec<u8>,
}

pub(crate) fn open_maybe_gzip_reader(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path).map_err(|e| AppError::io(format!("open {}", path.display()), e))?;
    if is_gzip_path(path)? {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

pub fn detect_input_format(path: &Path) -> Result<InputFormat> {
    if let Some(ext) = normalized_extension(path) {
        if ext == "fa" || ext == "fasta" || ext == "faa" || ext == "fna" {
            return Ok(InputFormat::Fasta);
        }
        if ext == "fq" || ext == "fastq" {
            return Ok(InputFormat::Fastq);
        }
    }

    let reader = open_maybe_gzip_reader(path)?;
    let mut reader = BufReader::new(reader);
    let mut buf = [0u8; 1];
    reader
        .read_exact(&mut buf)
        .map_err(|e| AppError::io(format!("read {}", path.display()), e))?;

    match buf[0] {
        b'>' => Ok(InputFormat::Fasta),
        b'@' => Ok(InputFormat::Fastq),
        other => Err(AppError::new(
            ErrorKind::Parse,
            format!("unknown input format: first byte {other}"),
        )),
    }
}

fn normalized_extension(path: &Path) -> Option<String> {
    let ext = path.extension()?.to_str()?.to_ascii_lowercase();
    if ext != "gz" {
        return Some(ext);
    }
    let stem = path.file_stem()?.to_str()?;
    let stem_path = Path::new(stem);
    stem_path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_ascii_lowercase())
}

fn is_gzip_path(path: &Path) -> Result<bool> {
    if path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        return Ok(true);
    }

    let mut f =
        File::open(path).map_err(|e| AppError::io(format!("open {}", path.display()), e))?;
    let mut magic = [0u8; 2];
    match f.read_exact(&mut magic) {
        Ok(_) => Ok(magic == [0x1f, 0x8b]),
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(false),
        Err(e) => Err(AppError::io(format!("read {}", path.display()), e)),
    }
}

pub fn validate_seq(seq: &[u8], dbtype: DbType) -> Result<()> {
    match dbtype {
        DbType::Auto => Ok(()),
        DbType::Protein => {
            if seq.iter().all(|b| b.is_ascii_uppercase()) {
                Ok(())
            } else {
                Err(AppError::new(
                    ErrorKind::Validation,
                    "protein sequence contains non A-Z characters",
                ))
            }
        }
        DbType::Nucleotide => {
            const ALLOWED: &[u8] = b"ACGTUNRYSWKMBDHV";
            if seq.iter().all(|b| ALLOWED.contains(b)) {
                Ok(())
            } else {
                Err(AppError::new(
                    ErrorKind::Validation,
                    "nucleotide sequence contains invalid characters",
                ))
            }
        }
    }
}
