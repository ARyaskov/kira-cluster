use std::collections::BTreeSet;
use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

fn run(mode: &str, base: &std::path::Path) -> String {
    let result = base.join(format!("result_{mode}"));
    let tmp = base.join(format!("tmp_{mode}"));

    let out = Command::new(bin())
        .args([
            "easy-cluster",
            "tests/data/cluster_sensitive.fasta",
            result.to_str().expect("utf8"),
            tmp.to_str().expect("utf8"),
            "--threads",
            "2",
            "--min-seq-id",
            "0.9",
            "--cov",
            "0.8",
            "--k-mer-size",
            "5",
            "--alignment-mode",
            mode,
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
    let mut set = BTreeSet::new();
    for line in tsv.lines() {
        if let Some(cid) = line.split('\t').next() {
            set.insert(cid.to_string());
        }
    }
    set.len()
}

#[test]
fn sensitive_mode_improves_recall_vs_fast() {
    let base = std::env::temp_dir().join("kira_cluster_stage3_sensitive");
    if base.exists() {
        std::fs::remove_dir_all(&base).expect("cleanup");
    }

    let fast = run("fast", &base);
    let sensitive = run("sensitive", &base);

    let fast_n = cluster_count(&fast);
    let sensitive_n = cluster_count(&sensitive);

    assert!(
        sensitive_n < fast_n,
        "fast:\n{fast}\nsensitive:\n{sensitive}"
    );
}
