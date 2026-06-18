use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

#[test]
fn compressed_and_legacy_segments_return_same_hits() {
    let base = std::env::temp_dir().join("kira_stage8_parity");
    let db = base.join("db");
    let idx_comp = base.join("idx_comp");
    let idx_legacy = base.join("idx_legacy");
    let out_comp = base.join("comp.tsv");
    let out_legacy = base.join("legacy.tsv");
    let _ = std::fs::remove_dir_all(&base);

    assert!(
        Command::new(bin())
            .args(["createdb", "tests/data/small.fasta", db.to_str().unwrap()])
            .output()
            .unwrap()
            .status
            .success()
    );
    assert!(
        Command::new(bin())
            .args([
                "build-index",
                db.to_str().unwrap(),
                idx_comp.to_str().unwrap(),
                "--k",
                "5",
                "--m",
                "8"
            ])
            .output()
            .unwrap()
            .status
            .success()
    );
    assert!(
        Command::new(bin())
            .args([
                "build-index",
                db.to_str().unwrap(),
                idx_legacy.to_str().unwrap(),
                "--k",
                "5",
                "--m",
                "8",
                "--legacy-uncompressed"
            ])
            .output()
            .unwrap()
            .status
            .success()
    );

    assert!(
        Command::new(bin())
            .args([
                "search-index",
                idx_comp.to_str().unwrap(),
                "tests/data/index_query.fasta",
                out_comp.to_str().unwrap(),
                "--k",
                "5",
                "--m",
                "8",
                "--top-k",
                "5"
            ])
            .output()
            .unwrap()
            .status
            .success()
    );
    assert!(
        Command::new(bin())
            .args([
                "search-index",
                idx_legacy.to_str().unwrap(),
                "tests/data/index_query.fasta",
                out_legacy.to_str().unwrap(),
                "--k",
                "5",
                "--m",
                "8",
                "--top-k",
                "5"
            ])
            .output()
            .unwrap()
            .status
            .success()
    );

    assert_eq!(
        std::fs::read(out_comp).unwrap(),
        std::fs::read(out_legacy).unwrap()
    );
}
