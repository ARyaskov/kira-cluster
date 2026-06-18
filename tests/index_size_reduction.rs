use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

fn dir_size(path: &std::path::Path) -> u64 {
    let mut total = 0u64;
    for e in walkdir(path) {
        if let Ok(m) = std::fs::metadata(&e) {
            if m.is_file() {
                total += m.len();
            }
        }
    }
    total
}

fn file_size(path: &std::path::Path) -> u64 {
    std::fs::metadata(path).map(|m| m.len()).unwrap_or(0)
}

fn walkdir(path: &std::path::Path) -> Vec<std::path::PathBuf> {
    let mut out = Vec::new();
    if let Ok(rd) = std::fs::read_dir(path) {
        for e in rd.flatten() {
            let p = e.path();
            if p.is_dir() {
                out.extend(walkdir(&p));
            } else {
                out.push(p);
            }
        }
    }
    out
}

#[test]
fn compressed_index_is_smaller_than_legacy() {
    let base = std::env::temp_dir().join("kira_stage8_size");
    let db = base.join("db");
    let idx_comp = base.join("idx_comp");
    let idx_legacy = base.join("idx_legacy");
    let _ = std::fs::remove_dir_all(&base);

    assert!(
        Command::new(bin())
            .args([
                "createdb",
                "tests/data/cluster_small.fasta",
                db.to_str().unwrap()
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
                idx_comp.to_str().unwrap(),
                "--k",
                "5",
                "--m",
                "12"
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
                "12",
                "--legacy-uncompressed"
            ])
            .output()
            .unwrap()
            .status
            .success()
    );

    let comp_seg = idx_comp.join("segment_0000");
    let leg_seg = idx_legacy.join("segment_0000");
    let comp_postings = file_size(&comp_seg.join("postings_hot.bin"))
        + file_size(&comp_seg.join("postings_cold.bin"));
    let legacy_postings = file_size(&leg_seg.join("postings.u32"));
    assert!(
        comp_postings < legacy_postings,
        "compressed_postings={} legacy_postings={}",
        comp_postings,
        legacy_postings
    );

    // End-to-end segment size may be larger on tiny fixtures due to extra Stage 8 metadata.
    // Keep a sanity bound to avoid accidental blow-ups.
    let s_comp = dir_size(&idx_comp);
    let s_leg = dir_size(&idx_legacy);
    assert!(
        s_comp <= s_leg + 1024,
        "compressed={} legacy={}",
        s_comp,
        s_leg
    );
}
