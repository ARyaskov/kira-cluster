use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

#[test]
fn profile_counters_are_sane() {
    let base = std::env::temp_dir().join("kira_profile_counter_sanity");
    let result = base.join("result");
    let tmp = base.join("tmp");
    let profile = base.join("profile.json");
    let _ = std::fs::remove_dir_all(&base);

    let out = Command::new(bin())
        .args([
            "easy-cluster",
            "tests/data/cluster_small.fasta",
            result.to_str().unwrap(),
            tmp.to_str().unwrap(),
            "--profile-json",
            profile.to_str().unwrap(),
        ])
        .output()
        .expect("run easy-cluster with profile");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let json = std::fs::read_to_string(&profile).expect("read profile json");
    let v: serde_json::Value = serde_json::from_str(&json).expect("parse profile json");

    assert!(v["n_sequences"].as_u64().unwrap_or(0) > 0);
    assert!(v["total_seeds_emitted"].as_u64().unwrap_or(0) > 0);
    assert!(v["unique_keys"].as_u64().unwrap_or(0) > 0);
    assert!(v["candidate_pairs_evaluated"].as_u64().unwrap_or(0) > 0);

    let cand = v["candidate_pairs_evaluated"].as_u64().unwrap_or(0);
    let h = v["hamming_rejected"].as_u64().unwrap_or(0);
    let u = v["ungapped_rejected"].as_u64().unwrap_or(0);
    let g = v["gapped_rejected"].as_u64().unwrap_or(0);
    assert!(h <= cand);
    assert!(u <= cand);
    assert!(g <= cand);

    assert!(v["stage_timings_ns"]["parse_ingest"].is_number());
    assert!(v["stage_timings_ns"]["seed_generation"].is_number());
    assert!(v["stage_timings_ns"]["final_write"].is_number());
}
