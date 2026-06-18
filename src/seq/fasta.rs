use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::error::{AppError, ErrorKind, Result};
use crate::seq::open_maybe_gzip_reader;

pub fn parse_fasta<F>(path: &Path, mut on_record: F) -> Result<()>
where
    F: FnMut(String, Vec<u8>) -> Result<()>,
{
    let reader = open_maybe_gzip_reader(path)?;
    let mut reader = BufReader::new(reader);
    let mut line = String::new();
    let mut current_name: Option<String> = None;
    let mut seq_buf: Vec<u8> = Vec::new();

    loop {
        line.clear();
        let n = reader
            .read_line(&mut line)
            .map_err(|e| AppError::io(format!("read {}", path.display()), e))?;
        if n == 0 {
            break;
        }

        let trimmed = line.trim_end_matches(['\n', '\r']);
        if trimmed.is_empty() {
            continue;
        }

        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(name) = current_name.take() {
                on_record(name, std::mem::take(&mut seq_buf))?;
            }
            current_name = Some(rest.trim().to_string());
            continue;
        }

        if current_name.is_none() {
            return Err(AppError::new(
                ErrorKind::Parse,
                "FASTA content before first header",
            ));
        }

        for &b in trimmed.as_bytes() {
            if !b.is_ascii_whitespace() {
                seq_buf.push(b.to_ascii_uppercase());
            }
        }
    }

    if let Some(name) = current_name {
        on_record(name, seq_buf)?;
    }

    Ok(())
}
