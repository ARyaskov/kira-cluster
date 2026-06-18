use std::path::Path;

use kira_fastq::{FastqFormat, FastqReader, ValidationMode};

use crate::error::{AppError, ErrorKind, Result};

pub fn parse_fastq<F>(path: &Path, mut on_record: F) -> Result<()>
where
    F: FnMut(String, Vec<u8>) -> Result<()>,
{
    let mut reader = FastqReader::from_path_auto(path)
        .map_err(|e| {
            AppError::new(
                ErrorKind::Parse,
                format!("parse FASTQ {}: {e}", path.display()),
            )
        })?
        .with_format(FastqFormat::MultiLine)
        .with_validation(ValidationMode::Qualities);

    while let Some(record) = reader.next().map_err(|e| {
        AppError::new(
            ErrorKind::Parse,
            format!("parse FASTQ {}: {e}", path.display()),
        )
    })? {
        let name = std::str::from_utf8(record.header())
            .map_err(|e| {
                AppError::new(ErrorKind::Parse, format!("FASTQ header is not UTF-8: {e}"))
            })?
            .trim()
            .to_string();
        let seq = record
            .seq()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect::<Vec<_>>();
        on_record(name, seq)?;
    }

    Ok(())
}
