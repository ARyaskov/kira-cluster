use crate::simd::SimdBackend;

#[derive(Debug, Clone, Copy)]
pub struct LogConfig {
    pub quiet: bool,
    pub verbose: bool,
}

impl LogConfig {
    pub fn new(quiet: bool, verbose: bool) -> Self {
        Self { quiet, verbose }
    }
}

pub fn info(cfg: LogConfig, event: &str, message: &str) {
    if !cfg.quiet {
        eprintln!("INFO {event} {message}");
    }
}

pub fn debug(cfg: LogConfig, event: &str, message: &str) {
    if cfg.verbose && !cfg.quiet {
        eprintln!("DEBUG {event} {message}");
    }
}

pub fn warn(cfg: LogConfig, code: &str, message: &str) {
    if !cfg.quiet {
        eprintln!("WARN {code} {message}");
    }
}

pub fn backend_startup(cfg: LogConfig, backend: SimdBackend) {
    let _ = cfg;
    eprintln!("INFO SIMD_BACKEND {}", backend.as_str());
}
