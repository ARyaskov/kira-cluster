use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

fn run(out_dir: &std::path::Path, threads: &str) {
    let result = out_dir.join("result");
    let tmp = out_dir.join("tmp");
    if out_dir.exists() {
        std::fs::remove_dir_all(out_dir).expect("cleanup");
    }

    let out = Command::new(bin())
        .args([
            "easy-cluster",
            "tests/data/cluster_small.fasta",
            result.to_str().expect("utf8"),
            tmp.to_str().expect("utf8"),
            "--threads",
            threads,
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
fn easy_cluster_is_deterministic_across_threads() {
    let base = std::env::temp_dir().join("kira_cluster_stage2_determinism");
    let a = base.join("a");
    let b = base.join("b");

    run(&a, "1");
    run(&b, "4");

    for name in ["clusters.tsv", "rep_seqs.fasta"] {
        let ba = std::fs::read(a.join("result").join(name)).expect("read a output");
        let bb = std::fs::read(b.join("result").join(name)).expect("read b output");
        assert_eq!(ba, bb, "mismatch in {name}");
    }
}
