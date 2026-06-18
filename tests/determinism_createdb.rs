use std::collections::BTreeSet;
use std::path::Path;
use std::process::Command;

fn bin() -> String {
    env!("CARGO_BIN_EXE_kira-cluster").to_string()
}

fn run_createdb(out_dir: &Path) {
    if out_dir.exists() {
        std::fs::remove_dir_all(out_dir).expect("cleanup");
    }

    let out = Command::new(bin())
        .args([
            "createdb",
            "tests/data/small.fasta",
            out_dir.to_str().expect("utf8 path"),
            "--gpu",
        ])
        .output()
        .expect("run createdb");

    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("WARN KW2001 GPU_FLAG_REQUESTED"));
}

#[test]
fn createdb_is_deterministic() {
    let base = std::env::temp_dir().join("kira_cluster_determinism");
    let a = base.join("a");
    let b = base.join("b");

    run_createdb(&a);
    run_createdb(&b);

    let names_a: BTreeSet<_> = std::fs::read_dir(&a)
        .expect("read a")
        .map(|e| e.expect("entry").file_name())
        .collect();
    let names_b: BTreeSet<_> = std::fs::read_dir(&b)
        .expect("read b")
        .map(|e| e.expect("entry").file_name())
        .collect();
    assert_eq!(names_a, names_b);

    for name in names_a {
        let pa = a.join(&name);
        let pb = b.join(&name);
        let ba = std::fs::read(&pa).expect("read pa");
        let bb = std::fs::read(&pb).expect("read pb");
        assert_eq!(ba, bb, "mismatch for {}", pa.display());
    }
}
