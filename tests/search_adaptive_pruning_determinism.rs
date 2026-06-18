use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

#[test]
fn adaptive_pruning_is_deterministic() {
    let base = std::env::temp_dir().join("kira_stage9_search_det");
    let db = base.join("db");
    let idx = base.join("idx");
    let out_a = base.join("out_a.tsv");
    let out_b = base.join("out_b.tsv");
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

    let run_search = |out: &std::path::Path| {
        Command::new(bin())
            .args([
                "search-index",
                idx.to_str().unwrap(),
                "tests/data/index_query.fasta",
                out.to_str().unwrap(),
                "--prune-df-quantile",
                "0.85",
                "--max-seeds-per-query",
                "64",
                "--work-budget",
                "10000",
            ])
            .status()
            .unwrap()
            .success()
    };

    assert!(run_search(&out_a));
    assert!(run_search(&out_b));

    let a = std::fs::read(&out_a).unwrap();
    let b = std::fs::read(&out_b).unwrap();
    assert_eq!(a, b);
}
