#[cfg(feature = "cuda")]
#[test]
fn gpu_matches_cpu_when_cuda_available() {
    use std::process::Command;

    fn bin() -> String {
        env!("CARGO_BIN_EXE_kira-cluster").to_string()
    }

    let base = std::env::temp_dir().join("kira_cluster_stage4_gpu_eq");
    if base.exists() {
        std::fs::remove_dir_all(&base).expect("cleanup");
    }

    let cpu_result = base.join("cpu_result");
    let cpu_tmp = base.join("cpu_tmp");
    let gpu_result = base.join("gpu_result");
    let gpu_tmp = base.join("gpu_tmp");

    let cpu = Command::new(bin())
        .args([
            "easy-cluster",
            "tests/data/cluster_small.fasta",
            cpu_result.to_str().expect("utf8"),
            cpu_tmp.to_str().expect("utf8"),
            "--cascade-level",
            "sensitive",
        ])
        .output()
        .expect("run cpu");
    assert!(
        cpu.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&cpu.stderr)
    );

    let gpu = Command::new(bin())
        .args([
            "easy-cluster",
            "tests/data/cluster_small.fasta",
            gpu_result.to_str().expect("utf8"),
            gpu_tmp.to_str().expect("utf8"),
            "--cascade-level",
            "sensitive",
            "--gpu",
            "--gpu-backend",
            "cuda",
        ])
        .output()
        .expect("run gpu");
    assert!(
        gpu.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&gpu.stderr)
    );

    let stderr = String::from_utf8_lossy(&gpu.stderr);
    if stderr.contains("GPU_BACKEND_UNAVAILABLE_FALLBACK") {
        return;
    }

    for name in ["clusters.tsv", "rep_seqs.fasta"] {
        let a = std::fs::read(cpu_result.join(name)).expect("read cpu");
        let b = std::fs::read(gpu_result.join(name)).expect("read gpu");
        assert_eq!(a, b, "mismatch in {name}");
    }
}
