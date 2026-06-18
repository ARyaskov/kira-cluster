use crate::error::{AppError, ErrorKind, Result};

pub fn not_available() -> Result<()> {
    Err(AppError::new(
        ErrorKind::Unsupported,
        "RMI_LITE_NOT_IMPLEMENTED",
    ))
}
