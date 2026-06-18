use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

fn run_profile(path: &std::path::Path, suffix: &str) {
    let base = std::env::temp_dir().join(format!("kira_profile_det_{suffix}"));
    let result = base.join("result");
    let tmp = base.join("tmp");
    let _ = std::fs::remove_dir_all(&base);

    let out = Command::new(bin())
        .args([
            "easy-cluster",
            "tests/data/cluster_small.fasta",
            result.to_str().unwrap(),
            tmp.to_str().unwrap(),
            "--profile-json",
            path.to_str().unwrap(),
        ])
        .output()
        .expect("run easy-cluster with profile-json");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn profile_json_is_deterministic_except_wallclock() {
    let base = std::env::temp_dir().join("kira_profile_det_files");
    let _ = std::fs::remove_dir_all(&base);
    std::fs::create_dir_all(&base).expect("create base");

    let a = base.join("a.json");
    let b = base.join("b.json");

    run_profile(&a, "a");
    run_profile(&b, "b");

    let mut va: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&a).expect("read a")).expect("parse a");
    let mut vb: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&b).expect("read b")).expect("parse b");

    // Wall-clock timing naturally varies across runs.
    va["total_wall_clock_ns"] = serde_json::Value::from(0u64);
    vb["total_wall_clock_ns"] = serde_json::Value::from(0u64);
    if let Some(obj) = va["stage_timings_ns"].as_object_mut() {
        for (_, v) in obj.iter_mut() {
            *v = serde_json::Value::from(0u64);
        }
    }
    if let Some(obj) = vb["stage_timings_ns"].as_object_mut() {
        for (_, v) in obj.iter_mut() {
            *v = serde_json::Value::from(0u64);
        }
    }

    assert_eq!(va, vb);
}
