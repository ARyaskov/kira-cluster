use std::fmt::{Display, Formatter};

pub type Result<T> = std::result::Result<T, AppError>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorKind {
    Cli,
    Io,
    Parse,
    Validation,
    Unsupported,
    Internal,
}

impl ErrorKind {
    pub fn code(self) -> &'static str {
        match self {
            ErrorKind::Cli => "KC1001",
            ErrorKind::Io => "KC1002",
            ErrorKind::Parse => "KC1003",
            ErrorKind::Validation => "KC1004",
            ErrorKind::Unsupported => "KC1005",
            ErrorKind::Internal => "KC1006",
        }
    }

    pub fn exit_status(self) -> i32 {
        match self {
            ErrorKind::Cli => 2,
            ErrorKind::Io => 3,
            ErrorKind::Parse => 4,
            ErrorKind::Validation => 5,
            ErrorKind::Unsupported => 6,
            ErrorKind::Internal => 10,
        }
    }
}

#[derive(Debug)]
pub struct AppError {
    pub kind: ErrorKind,
    pub message: String,
}

impl AppError {
    pub fn new(kind: ErrorKind, message: impl Into<String>) -> Self {
        Self {
            kind,
            message: message.into(),
        }
    }

    pub fn io(context: impl Into<String>, err: std::io::Error) -> Self {
        Self::new(ErrorKind::Io, format!("{}: {}", context.into(), err))
    }

    pub fn code(&self) -> &'static str {
        self.kind.code()
    }

    pub fn exit_status(&self) -> i32 {
        self.kind.exit_status()
    }
}

impl Display for AppError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {}", self.code(), self.message)
    }
}

impl std::error::Error for AppError {}

impl From<std::io::Error> for AppError {
    fn from(value: std::io::Error) -> Self {
        Self::io("io error", value)
    }
}
