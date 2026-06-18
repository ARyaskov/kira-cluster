use std::collections::BTreeSet;
use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

#[test]
fn easy_cluster_produces_expected_cluster_count() {
    let base = std::env::temp_dir().join("kira_cluster_stage2_small");
    let result = base.join("result");
    let tmp = base.join("tmp");
    if base.exists() {
        std::fs::remove_dir_all(&base).expect("cleanup");
    }

    let out = Command::new(bin())
        .args([
            "easy-cluster",
            "tests/data/cluster_small.fasta",
            result.to_str().expect("utf8"),
            tmp.to_str().expect("utf8"),
            "--threads",
            "2",
            "--min-seq-id",
            "0.9",
            "--cov",
            "0.8",
        ])
        .output()
        .expect("run easy-cluster");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let clusters = std::fs::read_to_string(result.join("clusters.tsv")).expect("read clusters.tsv");
    let mut ids = BTreeSet::new();
    let mut lines = 0usize;
    for line in clusters.lines() {
        let mut parts = line.split('\t');
        let cid = parts.next().expect("cid").parse::<u32>().expect("cid u32");
        ids.insert(cid);
        lines += 1;
    }

    assert_eq!(ids.len(), 2, "clusters.tsv:\n{clusters}");
    assert_eq!(lines, 5, "clusters.tsv:\n{clusters}");
    assert!(result.join("rep_seqs.fasta").exists());
}
