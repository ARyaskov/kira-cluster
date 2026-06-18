use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

#[test]
fn update_index_appends_manifest_segment() {
    let base = std::env::temp_dir().join("kira_stage7_update");
    let db1 = base.join("db1");
    let db2 = base.join("db2");
    let idx = base.join("idx");

    let _ = std::fs::remove_dir_all(&base);

    let c1 = Command::new(bin())
        .args([
            "createdb",
            "tests/data/small.fasta",
            db1.to_str().expect("utf8"),
        ])
        .output()
        .expect("run createdb 1");
    assert!(c1.status.success());

    let c2 = Command::new(bin())
        .args([
            "createdb",
            "tests/data/cluster_small.fasta",
            db2.to_str().expect("utf8"),
        ])
        .output()
        .expect("run createdb 2");
    assert!(c2.status.success());

    let b = Command::new(bin())
        .args([
            "build-index",
            db1.to_str().expect("utf8"),
            idx.to_str().expect("utf8"),
            "--k",
            "5",
            "--m",
            "8",
        ])
        .output()
        .expect("run build-index");
    assert!(b.status.success());

    let u = Command::new(bin())
        .args([
            "update-index",
            idx.to_str().expect("utf8"),
            db2.to_str().expect("utf8"),
            "--k",
            "5",
            "--m",
            "8",
        ])
        .output()
        .expect("run update-index");
    assert!(
        u.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&u.stderr)
    );

    let manifest = std::fs::read_to_string(idx.join("manifest.json")).expect("read manifest");
    assert!(manifest.contains("\"segment_id\": 0"));
    assert!(manifest.contains("\"segment_id\": 1"));

    let seg0 = idx.join("segment_0000");
    let seg1 = idx.join("segment_0001");
    assert!(seg0.join("kv_pgm.bin").exists());
    assert!(seg1.join("kv_pgm.bin").exists());

    let meta0: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(seg0.join("segment_meta.json")).unwrap())
            .unwrap();
    let meta1: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(seg1.join("segment_meta.json")).unwrap())
            .unwrap();
    assert_eq!(meta0["global_seq_base"].as_u64().unwrap(), 0);
    assert_eq!(
        meta1["global_seq_base"].as_u64().unwrap(),
        meta0["n_seqs"].as_u64().unwrap()
    );

    let query = base.join("segment1_query.fasta");
    let out_tsv = base.join("segment1_hits.tsv");
    std::fs::write(&query, b">q_segment1\nMKTAYIAKQRQISFVKSHFSRQDILDLWQ\n").unwrap();

    let s = Command::new(bin())
        .args([
            "search-index",
            idx.to_str().expect("utf8"),
            query.to_str().expect("utf8"),
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
    assert!(
        s.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&s.stderr)
    );

    let hits = std::fs::read_to_string(out_tsv).expect("read hits");
    let first = hits.lines().next().expect("first hit");
    let cols: Vec<_> = first.split('\t').collect();
    assert_eq!(cols[0], "q_segment1");
    assert_eq!(cols[1], meta0["n_seqs"].as_u64().unwrap().to_string());
    assert_eq!(cols[6], "1");
    assert_eq!(cols[8], "seq1");
    assert_eq!(cols[9], "0");
}
