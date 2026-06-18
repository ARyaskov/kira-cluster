use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

#[test]
fn pgm_files_are_deterministic() {
    let base = std::env::temp_dir().join("kira_stage9_pgm_det");
    let db = base.join("db");
    let idx_a = base.join("idx_a");
    let idx_b = base.join("idx_b");
    let _ = std::fs::remove_dir_all(&base);

    assert!(
        Command::new(bin())
            .args(["createdb", "tests/data/small.fasta", db.to_str().unwrap()])
            .status()
            .unwrap()
            .success()
    );

    for out in [&idx_a, &idx_b] {
        assert!(
            Command::new(bin())
                .args([
                    "build-index",
                    db.to_str().unwrap(),
                    out.to_str().unwrap(),
                    "--k",
                    "5",
                    "--m",
                    "8",
                    "--pgm-epsilon",
                    "64",
                ])
                .status()
                .unwrap()
                .success()
        );
    }

    for name in ["segment_0000/kv_pgm.bin"] {
        let a = std::fs::read(idx_a.join(name)).unwrap();
        let b = std::fs::read(idx_b.join(name)).unwrap();
        assert_eq!(a, b, "mismatch for {name}");
    }
}
