use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

#[test]
fn build_index_is_byte_deterministic() {
    let base = std::env::temp_dir().join("kira_stage8_det");
    let db = base.join("db");
    let a = base.join("a");
    let b = base.join("b");
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
                a.to_str().unwrap(),
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
                b.to_str().unwrap(),
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

    for f in [
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
    ] {
        assert_eq!(
            std::fs::read(a.join(f)).unwrap(),
            std::fs::read(b.join(f)).unwrap(),
            "{f}"
        );
    }
}
