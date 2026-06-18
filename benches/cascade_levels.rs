use std::process::Command;
use std::time::Instant;

fn kira_bin() -> String {
    std::env::var("KIRA_CLUSTER_BIN").unwrap_or_else(|_| "target/release/kira-cluster".to_string())
}

fn run(level: &str) -> std::time::Duration {
    let base = std::env::temp_dir().join(format!("kira_bench_{level}"));
    let result = base.join("result");
    let tmp = base.join("tmp");
    let _ = std::fs::remove_dir_all(&base);

    let t0 = Instant::now();
    let out = Command::new(kira_bin())
        .args([
            "easy-cluster",
            "tests/data/cluster_small.fasta",
            result.to_str().expect("utf8"),
            tmp.to_str().expect("utf8"),
            "--cascade-level",
            level,
            "--batch-size",
            "64",
        ])
        .output()
        .expect("run easy-cluster");

    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    t0.elapsed()
}

fn main() {
    let levels = ["fast", "medium", "sensitive", "ultra"];
    println!("level\twall_ms");
    for level in levels {
        let elapsed = run(level);
        println!("{level}\t{}", elapsed.as_secs_f64() * 1000.0);
    }
}
