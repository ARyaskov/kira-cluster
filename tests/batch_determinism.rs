use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

fn run(base: &std::path::Path, batch_size: &str) {
    let result = base.join(format!("result_{batch_size}"));
    let tmp = base.join(format!("tmp_{batch_size}"));

    let out = Command::new(bin())
        .args([
            "easy-cluster",
            "tests/data/cluster_small.fasta",
            result.to_str().expect("utf8"),
            tmp.to_str().expect("utf8"),
            "--cascade-level",
            "sensitive",
            "--batch-size",
            batch_size,
            "--threads",
            "4",
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
}

#[test]
fn outputs_are_identical_across_batch_sizes() {
    let base = std::env::temp_dir().join("kira_cluster_stage4_batch_det");
    if base.exists() {
        std::fs::remove_dir_all(&base).expect("cleanup");
    }

    run(&base, "3");
    run(&base, "64");

    for name in ["clusters.tsv", "rep_seqs.fasta"] {
        let a = std::fs::read(base.join("result_3").join(name)).expect("read a");
        let b = std::fs::read(base.join("result_64").join(name)).expect("read b");
        assert_eq!(a, b, "mismatch in {name}");
    }
}
