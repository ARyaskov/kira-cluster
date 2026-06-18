use std::collections::BTreeSet;
use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

fn run(level: &str, base: &std::path::Path) -> String {
    let result = base.join(format!("result_{level}"));
    let tmp = base.join(format!("tmp_{level}"));

    let out = Command::new(bin())
        .args([
            "easy-cluster",
            "tests/data/cluster_sensitive.fasta",
            result.to_str().expect("utf8"),
            tmp.to_str().expect("utf8"),
            "--cascade-level",
            level,
            "--min-seq-id",
            "0.9",
            "--cov",
            "0.8",
            "--k-mer-size",
            "5",
        ])
        .output()
        .expect("run easy-cluster");

    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    std::fs::read_to_string(result.join("clusters.tsv")).expect("read clusters.tsv")
}

fn cluster_count(tsv: &str) -> usize {
    let mut ids = BTreeSet::new();
    for line in tsv.lines() {
        if let Some(cid) = line.split('\t').next() {
            ids.insert(cid.to_string());
        }
    }
    ids.len()
}

#[test]
fn cascade_ultra_is_not_worse_than_fast() {
    let base = std::env::temp_dir().join("kira_cluster_stage4_cascade");
    if base.exists() {
        std::fs::remove_dir_all(&base).expect("cleanup");
    }

    let fast = run("fast", &base);
    let ultra = run("ultra", &base);
    assert!(cluster_count(&ultra) <= cluster_count(&fast));
}
