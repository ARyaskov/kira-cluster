mod alignment;
mod cascade;
mod cli;
mod cluster;
mod db;
mod error;
mod gpu;
mod index;
mod io;
mod profile;
mod scheduler;
mod seq;
mod simd;
mod util;

fn main() {
    if let Err(err) = cli::run() {
        eprintln!("ERROR {}", err);
        std::process::exit(err.exit_status());
    }
}
