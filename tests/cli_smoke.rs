use std::process::Command;
use std::{fs::File, io::Write};

use flate2::Compression;
use flate2::write::GzEncoder;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

#[test]
fn help_works() {
    let out = Command::new(bin())
        .arg("--help")
        .output()
        .expect("run --help");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn createdb_works() {
    let root = std::env::temp_dir().join("kira_cluster_cli_smoke_db");
    if root.exists() {
        std::fs::remove_dir_all(&root).expect("cleanup root");
    }

    let out = Command::new(bin())
        .args([
            "createdb",
            "tests/data/small.fasta",
            root.to_str().expect("utf8 path"),
        ])
        .output()
        .expect("run createdb");

    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(root.join("meta.json").exists());
    assert!(root.join("seqs.bin").exists());
}

#[test]
fn createdb_fastq_gz_works() {
    let base = std::env::temp_dir().join("kira_cluster_cli_smoke_fastq_gz");
    let input = base.join("small.fastq.gz");
    let out_db = base.join("db");
    if base.exists() {
        std::fs::remove_dir_all(&base).expect("cleanup base");
    }
    std::fs::create_dir_all(&base).expect("create base");

    let raw_fastq = std::fs::read("tests/data/small.fastq").expect("read small.fastq");
    let mut enc = GzEncoder::new(
        File::create(&input).expect("create fastq.gz"),
        Compression::default(),
    );
    enc.write_all(&raw_fastq).expect("write gz payload");
    enc.finish().expect("finish gz");

    let out = Command::new(bin())
        .args([
            "createdb",
            input.to_str().expect("utf8 path"),
            out_db.to_str().expect("utf8 path"),
        ])
        .output()
        .expect("run createdb fastq.gz");

    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(out_db.join("meta.json").exists());
    assert!(out_db.join("seqs.bin").exists());
}

#[test]
fn fastq_dbtype_controls_alphabet_validation() {
    let base = std::env::temp_dir().join("kira_cluster_cli_smoke_fastq_dbtype");
    let protein_db = base.join("protein_db");
    let nucleotide_db = base.join("nucleotide_db");
    if base.exists() {
        std::fs::remove_dir_all(&base).expect("cleanup base");
    }
    std::fs::create_dir_all(&base).expect("create base");

    let protein = Command::new(bin())
        .args([
            "createdb",
            "tests/data/small.fastq",
            protein_db.to_str().expect("utf8 path"),
            "--dbtype",
            "protein",
        ])
        .output()
        .expect("run protein createdb");
    assert!(
        protein.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&protein.stderr)
    );

    let nucleotide = Command::new(bin())
        .args([
            "createdb",
            "tests/data/small.fastq",
            nucleotide_db.to_str().expect("utf8 path"),
            "--dbtype",
            "nucleotide",
        ])
        .output()
        .expect("run nucleotide createdb");
    assert!(!nucleotide.status.success());
    assert!(
        String::from_utf8_lossy(&nucleotide.stderr)
            .contains("nucleotide sequence contains invalid characters")
    );
}
