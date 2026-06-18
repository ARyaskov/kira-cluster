use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

#[test]
fn index_search_is_deterministic() {
    let base = std::env::temp_dir().join("kira_stage7_search_det");
    let dbdir = base.join("db");
    let idx = base.join("idx");
    let out_a = base.join("a.tsv");
    let out_b = base.join("b.tsv");

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

    for out_tsv in [&out_a, &out_b] {
        let out = Command::new(bin())
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
                "5",
            ])
            .output()
            .expect("run search-index");
        assert!(
            out.status.success(),
            "stderr: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }

    let a = std::fs::read(&out_a).expect("read a");
    let b = std::fs::read(&out_b).expect("read b");
    assert_eq!(a, b);
}
