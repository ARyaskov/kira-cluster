use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsString;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

#[derive(Debug)]
struct BenchResult {
    tool: &'static str,
    status: &'static str,
    wall_ms: Option<f64>,
    peak_rss_kb: Option<u64>,
    clusters: Option<usize>,
    pairwise_recall_vs_kira: Option<f64>,
}

fn main() {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let input = args
        .first()
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("tests/data/cluster_small.fasta"));
    let out_dir = args
        .get(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("target/bench-tools"));

    fs::create_dir_all(&out_dir).expect("create bench output dir");
    let names = read_fasta_names(&input).expect("read input names");

    let kira = run_kira(&input, &out_dir, &names);
    let kira_clusters =
        load_kira_clusters(&out_dir.join("kira").join("clusters.tsv"), &names).unwrap_or_default();

    let mut results = vec![kira];
    results.push(run_mmseqs(&input, &out_dir, &kira_clusters));
    results.push(run_cdhit(&input, &out_dir, &kira_clusters));
    results.push(run_vsearch(&input, &out_dir, &kira_clusters));

    println!("tool\tstatus\twall_ms\tpeak_rss_kb\tclusters\tpairwise_recall_vs_kira");
    for r in results {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            r.tool,
            r.status,
            r.wall_ms
                .map(|v| format!("{v:.3}"))
                .unwrap_or_else(|| "NA".to_string()),
            r.peak_rss_kb
                .map(|v| v.to_string())
                .unwrap_or_else(|| "NA".to_string()),
            r.clusters
                .map(|v| v.to_string())
                .unwrap_or_else(|| "NA".to_string()),
            r.pairwise_recall_vs_kira
                .map(|v| format!("{v:.6}"))
                .unwrap_or_else(|| "NA".to_string())
        );
    }
}

fn kira_bin() -> OsString {
    std::env::var_os("KIRA_CLUSTER_BIN")
        .unwrap_or_else(|| OsString::from("target/release/kira-cluster"))
}

fn run_kira(input: &Path, out_dir: &Path, names: &[String]) -> BenchResult {
    let result = out_dir.join("kira");
    let tmp = out_dir.join("kira_tmp");
    let _ = fs::remove_dir_all(&result);
    let _ = fs::remove_dir_all(&tmp);

    let program = kira_bin();
    let args = vec![
        OsString::from("easy-cluster"),
        input.as_os_str().to_os_string(),
        result.as_os_str().to_os_string(),
        tmp.as_os_str().to_os_string(),
        OsString::from("--min-seq-id"),
        OsString::from("0.9"),
        OsString::from("--cov"),
        OsString::from("0.8"),
    ];
    let measured = run_measured(program, args);
    let clusters = load_kira_clusters(&result.join("clusters.tsv"), names).unwrap_or_default();
    BenchResult {
        tool: "kira-cluster",
        status: if measured.ok { "ok" } else { "failed" },
        wall_ms: Some(measured.wall_ms),
        peak_rss_kb: measured.peak_rss_kb,
        clusters: Some(clusters.len()),
        pairwise_recall_vs_kira: Some(1.0),
    }
}

fn run_mmseqs(input: &Path, out_dir: &Path, kira_clusters: &[Vec<String>]) -> BenchResult {
    if !exists_on_path("mmseqs") {
        return skipped("mmseqs2");
    }
    let prefix = out_dir.join("mmseqs");
    let tmp = out_dir.join("mmseqs_tmp");
    let cluster_tsv = out_dir.join("mmseqs_cluster.tsv");
    let _ = fs::remove_dir_all(&tmp);
    let _ = fs::remove_file(&cluster_tsv);

    let args = vec![
        OsString::from("easy-cluster"),
        input.as_os_str().to_os_string(),
        prefix.as_os_str().to_os_string(),
        tmp.as_os_str().to_os_string(),
        OsString::from("--min-seq-id"),
        OsString::from("0.9"),
        OsString::from("-c"),
        OsString::from("0.8"),
        OsString::from("--cov-mode"),
        OsString::from("1"),
    ];
    let measured = run_measured(OsString::from("mmseqs"), args);
    let clusters = load_mmseqs_clusters(&cluster_tsv).unwrap_or_default();
    measured_result("mmseqs2", measured, clusters, kira_clusters)
}

fn run_cdhit(input: &Path, out_dir: &Path, kira_clusters: &[Vec<String>]) -> BenchResult {
    if !exists_on_path("cd-hit") {
        return skipped("cd-hit");
    }
    let output = out_dir.join("cdhit");
    let _ = fs::remove_file(&output);
    let _ = fs::remove_file(output.with_extension("clstr"));
    let args = vec![
        OsString::from("-i"),
        input.as_os_str().to_os_string(),
        OsString::from("-o"),
        output.as_os_str().to_os_string(),
        OsString::from("-c"),
        OsString::from("0.9"),
        OsString::from("-n"),
        OsString::from("5"),
        OsString::from("-M"),
        OsString::from("0"),
        OsString::from("-T"),
        OsString::from("0"),
    ];
    let measured = run_measured(OsString::from("cd-hit"), args);
    let clusters = load_cdhit_clusters(&output.with_extension("clstr")).unwrap_or_default();
    measured_result("cd-hit", measured, clusters, kira_clusters)
}

fn run_vsearch(input: &Path, out_dir: &Path, kira_clusters: &[Vec<String>]) -> BenchResult {
    if !exists_on_path("vsearch") {
        return skipped("vsearch");
    }
    let centroids = out_dir.join("vsearch_centroids.fasta");
    let uc = out_dir.join("vsearch.uc");
    let _ = fs::remove_file(&centroids);
    let _ = fs::remove_file(&uc);
    let args = vec![
        OsString::from("--cluster_fast"),
        input.as_os_str().to_os_string(),
        OsString::from("--id"),
        OsString::from("0.9"),
        OsString::from("--centroids"),
        centroids.as_os_str().to_os_string(),
        OsString::from("--uc"),
        uc.as_os_str().to_os_string(),
    ];
    let measured = run_measured(OsString::from("vsearch"), args);
    let clusters = load_vsearch_clusters(&uc).unwrap_or_default();
    measured_result("vsearch", measured, clusters, kira_clusters)
}

fn measured_result(
    tool: &'static str,
    measured: Measured,
    clusters: Vec<Vec<String>>,
    kira_clusters: &[Vec<String>],
) -> BenchResult {
    BenchResult {
        tool,
        status: if measured.ok { "ok" } else { "failed" },
        wall_ms: Some(measured.wall_ms),
        peak_rss_kb: measured.peak_rss_kb,
        clusters: Some(clusters.len()),
        pairwise_recall_vs_kira: Some(pairwise_recall(kira_clusters, &clusters)),
    }
}

fn skipped(tool: &'static str) -> BenchResult {
    BenchResult {
        tool,
        status: "skipped",
        wall_ms: None,
        peak_rss_kb: None,
        clusters: None,
        pairwise_recall_vs_kira: None,
    }
}

struct Measured {
    ok: bool,
    wall_ms: f64,
    peak_rss_kb: Option<u64>,
}

fn run_measured(program: OsString, args: Vec<OsString>) -> Measured {
    if Path::new("/usr/bin/time").exists() {
        let mut time_args = Vec::new();
        if cfg!(target_os = "macos") {
            time_args.push(OsString::from("-l"));
        } else {
            time_args.push(OsString::from("-v"));
        }
        time_args.push(program.clone());
        time_args.extend(args.clone());

        let t0 = Instant::now();
        let output = Command::new("/usr/bin/time").args(time_args).output();
        if let Ok(output) = output {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Measured {
                ok: output.status.success(),
                wall_ms: t0.elapsed().as_secs_f64() * 1000.0,
                peak_rss_kb: parse_peak_rss_kb(&stderr),
            };
        }
    }

    let t0 = Instant::now();
    let output = Command::new(program).args(args).output();
    Measured {
        ok: output.map(|o| o.status.success()).unwrap_or(false),
        wall_ms: t0.elapsed().as_secs_f64() * 1000.0,
        peak_rss_kb: None,
    }
}

fn parse_peak_rss_kb(stderr: &str) -> Option<u64> {
    for line in stderr.lines() {
        let lower = line.to_ascii_lowercase();
        if lower.contains("maximum resident set size") {
            let n = line
                .split_whitespace()
                .find_map(|part| part.parse::<u64>().ok())?;
            return if cfg!(target_os = "macos") {
                Some(n / 1024)
            } else {
                Some(n)
            };
        }
    }
    None
}

fn exists_on_path(program: &str) -> bool {
    let Some(paths) = std::env::var_os("PATH") else {
        return false;
    };
    std::env::split_paths(&paths).any(|dir| dir.join(program).is_file())
}

fn read_fasta_names(path: &Path) -> std::io::Result<Vec<String>> {
    let file = fs::File::open(path)?;
    let reader = BufReader::new(file);
    let mut names = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix('>') {
            names.push(rest.split_whitespace().next().unwrap_or(rest).to_string());
        }
    }
    Ok(names)
}

fn load_kira_clusters(path: &Path, names: &[String]) -> std::io::Result<Vec<Vec<String>>> {
    let file = fs::File::open(path)?;
    let reader = BufReader::new(file);
    let mut groups: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for line in reader.lines() {
        let line = line?;
        let cols = line.split('\t').collect::<Vec<_>>();
        if cols.len() < 3 {
            continue;
        }
        let Some(member_id) = cols[2].parse::<usize>().ok() else {
            continue;
        };
        let Some(name) = names.get(member_id) else {
            continue;
        };
        groups
            .entry(cols[0].to_string())
            .or_default()
            .push(name.clone());
    }
    Ok(groups.into_values().collect())
}

fn load_mmseqs_clusters(path: &Path) -> std::io::Result<Vec<Vec<String>>> {
    let file = fs::File::open(path)?;
    let reader = BufReader::new(file);
    let mut groups: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for line in reader.lines() {
        let line = line?;
        let cols = line.split('\t').collect::<Vec<_>>();
        if cols.len() < 2 {
            continue;
        }
        groups
            .entry(cols[0].to_string())
            .or_default()
            .push(cols[1].to_string());
    }
    Ok(groups.into_values().collect())
}

fn load_cdhit_clusters(path: &Path) -> std::io::Result<Vec<Vec<String>>> {
    let file = fs::File::open(path)?;
    let reader = BufReader::new(file);
    let mut clusters = Vec::new();
    let mut current = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with(">Cluster") {
            if !current.is_empty() {
                clusters.push(std::mem::take(&mut current));
            }
            continue;
        }
        if let Some(start) = line.find('>') {
            let rest = &line[start + 1..];
            let end = rest.find("...").unwrap_or(rest.len());
            current.push(
                rest[..end]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string(),
            );
        }
    }
    if !current.is_empty() {
        clusters.push(current);
    }
    Ok(clusters)
}

fn load_vsearch_clusters(path: &Path) -> std::io::Result<Vec<Vec<String>>> {
    let file = fs::File::open(path)?;
    let reader = BufReader::new(file);
    let mut groups: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for line in reader.lines() {
        let line = line?;
        let cols = line.split('\t').collect::<Vec<_>>();
        if cols.len() < 9 {
            continue;
        }
        let rec_type = cols[0];
        if rec_type == "S" || rec_type == "H" {
            groups
                .entry(cols[1].to_string())
                .or_default()
                .push(cols[8].to_string());
        }
    }
    Ok(groups.into_values().collect())
}

fn pairwise_recall(reference: &[Vec<String>], candidate: &[Vec<String>]) -> f64 {
    let reference_pairs = cluster_pairs(reference);
    if reference_pairs.is_empty() {
        return 1.0;
    }
    let candidate_pairs = cluster_pairs(candidate);
    let hits = reference_pairs
        .iter()
        .filter(|pair| candidate_pairs.contains(*pair))
        .count();
    hits as f64 / reference_pairs.len() as f64
}

fn cluster_pairs(clusters: &[Vec<String>]) -> BTreeSet<(String, String)> {
    let mut out = BTreeSet::new();
    for cluster in clusters {
        for i in 0..cluster.len() {
            for j in i + 1..cluster.len() {
                let a = cluster[i].clone();
                let b = cluster[j].clone();
                if a <= b {
                    out.insert((a, b));
                } else {
                    out.insert((b, a));
                }
            }
        }
    }
    out
}
