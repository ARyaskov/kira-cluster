use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

fn run_createdb(dbdir: &std::path::Path) {
    if dbdir.exists() {
        std::fs::remove_dir_all(dbdir).expect("cleanup dbdir");
    }
    let out = Command::new(bin())
        .args([
            "createdb",
            "tests/data/small.fasta",
            dbdir.to_str().expect("utf8"),
        ])
        .output()
        .expect("run createdb");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

fn run_build_index(dbdir: &std::path::Path, indexdir: &std::path::Path) {
    if indexdir.exists() {
        std::fs::remove_dir_all(indexdir).expect("cleanup indexdir");
    }
    let out = Command::new(bin())
        .args([
            "build-index",
            dbdir.to_str().expect("utf8"),
            indexdir.to_str().expect("utf8"),
            "--k",
            "5",
            "--m",
            "8",
        ])
        .output()
        .expect("run build-index");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn index_build_is_deterministic() {
    let base = std::env::temp_dir().join("kira_stage7_build_det");
    let dbdir = base.join("db");
    let idx_a = base.join("idx_a");
    let idx_b = base.join("idx_b");

    run_createdb(&dbdir);
    run_build_index(&dbdir, &idx_a);
    run_build_index(&dbdir, &idx_b);

    for name in [
        "manifest.json",
        "segment_0000/segment_meta.json",
        "segment_0000/keys.u64",
        "segment_0000/df.u32",
        "segment_0000/tiers.u8",
        "segment_0000/tier_pos.u32",
        "segment_0000/offsets_hot.u64",
        "segment_0000/offsets_cold.u64",
        "segment_0000/postings_hot.bin",
        "segment_0000/postings_cold.bin",
        "segment_0000/codec_hot.json",
        "segment_0000/codec_cold.json",
        "segment_0000/bloom.bin",
        "segment_0000/bloom_meta.json",
        "segment_0000/kv_pgm.bin",
        "segment_0000/segment_stats.json",
        "seed_policy.json",
    ] {
        let a = std::fs::read(idx_a.join(name)).expect("read a");
        let b = std::fs::read(idx_b.join(name)).expect("read b");
        assert_eq!(a, b, "mismatch for {name}");
    }
}
