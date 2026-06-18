use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

#[test]
fn index_search_returns_expected_top_hit() {
    let base = std::env::temp_dir().join("kira_stage7_correctness");
    let dbdir = base.join("db");
    let idx = base.join("idx");
    let out_tsv = base.join("out.tsv");

    let _ = std::fs::remove_dir_all(&base);

    let createdb = Command::new(bin())
        .args([
            "createdb",
            "tests/data/small.fasta",
            dbdir.to_str().expect("utf8"),
        ])
        .output()
        .expect("run createdb");
    assert!(createdb.status.success());

    let build = Command::new(bin())
        .args([
            "build-index",
            dbdir.to_str().expect("utf8"),
            idx.to_str().expect("utf8"),
            "--k",
            "5",
            "--m",
            "8",
        ])
        .output()
        .expect("run build-index");
    assert!(build.status.success());

    let search = Command::new(bin())
        .args([
            "search-index",
            idx.to_str().expect("utf8"),
            "tests/data/index_query.fasta",
            out_tsv.to_str().expect("utf8"),
            "--k",
            "5",
            "--m",
            "8",
            "--top-k",
            "1",
        ])
        .output()
        .expect("run search-index");
    assert!(search.status.success());

    let tsv = std::fs::read_to_string(out_tsv).expect("read out");
    let first = tsv.lines().next().expect("first line");
    let cols: Vec<_> = first.split('\t').collect();
    assert_eq!(cols[0], "q1");
    assert_eq!(cols[1], "0");
}
