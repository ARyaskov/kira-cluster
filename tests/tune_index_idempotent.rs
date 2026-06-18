use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

#[test]
fn tune_index_is_idempotent() {
    let base = std::env::temp_dir().join("kira_stage9_tune_idempotent");
    let db = base.join("db");
    let idx = base.join("idx");
    let _ = std::fs::remove_dir_all(&base);

    assert!(
        Command::new(bin())
            .args(["createdb", "tests/data/small.fasta", db.to_str().unwrap()])
            .status()
            .unwrap()
            .success()
    );

    assert!(
        Command::new(bin())
            .args(["build-index", db.to_str().unwrap(), idx.to_str().unwrap()])
            .status()
            .unwrap()
            .success()
    );

    let run_tune = || {
        Command::new(bin())
            .args([
                "tune-index",
                idx.to_str().unwrap(),
                "--pgm-epsilon",
                "64",
                "--seed-policy",
                "variable-k",
            ])
            .status()
            .unwrap()
            .success()
    };

    assert!(run_tune());
    let a_meta = std::fs::read(idx.join("tune_meta.json")).unwrap();
    let a_seed = std::fs::read(idx.join("seed_policy.json")).unwrap();
    let a_pgm = std::fs::read(idx.join("segment_0000/kv_pgm.bin")).unwrap();

    assert!(run_tune());
    let b_meta = std::fs::read(idx.join("tune_meta.json")).unwrap();
    let b_seed = std::fs::read(idx.join("seed_policy.json")).unwrap();
    let b_pgm = std::fs::read(idx.join("segment_0000/kv_pgm.bin")).unwrap();

    assert_eq!(a_meta, b_meta);
    assert_eq!(a_seed, b_seed);
    assert_eq!(a_pgm, b_pgm);
}
