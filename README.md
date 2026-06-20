# kira-cluster

`kira-cluster` is a deterministic Rust CLI for MMseqs2-like approximate high-throughput sequence clustering/search workflows.

## Contract

The project contract is MMseqs2-like approximate high-throughput clustering:

- optimize for large FASTA/FASTQ inputs, deterministic output, and bounded work per query;
- use k-mer/minimizer candidate generation plus cascaded filters rather than exhaustive all-vs-all alignment;
- prefer stable, reproducible cluster representatives over exact transitive closure;
- expose MMseqs2-shaped commands where practical, but do not promise byte-for-byte MMseqs2 compatibility.

## Competitive Goal

`kira-cluster` is intentionally positioned as a small deterministic MMseqs2-like engine:

- reproducible clustering/search output from the same inputs, seed, and options;
- incremental multi-segment indexes with stable global sequence IDs;
- Rust-native embeddability and simple filesystem artifacts;
- high-throughput approximate clustering/search, not a full replacement for the mature MMseqs2/CD-HIT/VSEARCH feature surface.

## FASTQ Contract

- FASTQ framing, gzip/BGZF detection, multi-line records, and quality validation are delegated to `kira-fastq`.
- Sequence alphabet validation is owned by `createdb --dbtype`.
- `--dbtype nucleotide` accepts nucleotide/IUPAC symbols only.
- `--dbtype protein` accepts uppercase protein-like ASCII letters.
- `--dbtype auto` keeps ingestion permissive and resolves the DB type to nucleotide only when all records look nucleotide-like; otherwise it resolves to protein.

## Support

- `createdb <input> <dbdir>`
- `easy-cluster <input> <result_dir> <tmp_dir>` (deterministic greedy clustering)
- `build-index <dbdir> <indexdir>`
- `search-index <indexdir> <query_fasta> <out_tsv>`
- `update-index <indexdir> <new_dbdir>`
- `tune-index <indexdir>`
- `serve <indexdir>` (HTTP stub + core API)
- `kira_kv_engine` PGM key lookup (`kv_pgm.bin`) per segment
- Adaptive query-time pruning (`--prune-df-quantile`, `--max-seeds-per-query`, `--work-budget`)
- Seed policy metadata (`seed_policy.json`) with deterministic variable-k mode
- Adaptive cascade levels: `--cascade-level <fast|medium|sensitive|ultra>`
- Batch scheduler: `--batch-size`, `--gpu-memory-limit`, `--cpu-threads`
- Sensitive alignment controls: `--alignment-mode`, `--sensitivity`, `--gap-open`, `--gap-extend`
- GPU path (feature-gated): `--gpu --gpu-backend cuda`
- Runtime SIMD backend report (`avx2`, `neon`, `scalar`)
- Anti-diagonal SIMD Gotoh score kernel for AVX2/NEON with scalar tie tracking for deterministic `matches/aligned_len`
- FASTA parsing and FASTQ parsing through `kira-fastq` (plain, gzip/BGZF, multi-line FASTQ)
- Global sequence IDs across appended index segments (`global_seq_base` in `segment_meta.json`)
- Deterministic DB directory output
- Cluster outputs:
  - `clusters.tsv`
  - `rep_seqs.fasta`
- Search output TSV columns:
  `query_id`, `target_global_id`, `score`, `identity`, `coverage`, `seed_hits`, `segment_id`, `explain`, `target_name`, `target_local_id`

## Examples

```bash
cargo run -- createdb tests/data/small.fasta target/tmp/db
cargo run -- easy-cluster tests/data/small.fasta target/tmp/result target/tmp/tmp --print-plan
cargo run -- easy-cluster tests/data/cluster_small.fasta target/tmp/result target/tmp/tmp
cargo run -- easy-cluster tests/data/cluster_sensitive.fasta target/tmp/result target/tmp/tmp --alignment-mode sensitive
cargo run -- easy-cluster tests/data/cluster_small.fasta target/tmp/result target/tmp/tmp --cascade-level ultra --batch-size 2048
cargo run -- build-index target/tmp/db target/tmp/index --k 5 --m 8
cargo run -- tune-index target/tmp/index --pgm-epsilon 64 --seed-policy variable-k
cargo run -- search-index target/tmp/index tests/data/index_query.fasta target/tmp/hits.tsv --k 5 --m 8 --top-k 10 --prune-df-quantile 0.9 --max-seeds-per-query 128
cargo run -- update-index target/tmp/index target/tmp/db2 --k 5 --m 8
cargo bench --bench tool_compare -- tests/data/cluster_small.fasta target/bench-tools
```

## UniRef50 Trial

UniRef50 is large. Use a release build, a fast local SSD, and enough free disk for the compressed FASTA plus the uncompressed DB/index artifacts.

```bash
mkdir -p data target/uniref50
curl -L -o data/uniref50.fasta.gz \
  https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz

cargo build --release

./target/release/kira-cluster createdb \
  data/uniref50.fasta.gz target/uniref50/db \
  --dbtype protein --threads 8

./target/release/kira-cluster easy-cluster \
  target/uniref50/db target/uniref50/clusters target/uniref50/tmp \
  --dbtype protein --min-seq-id 0.50 --cov 0.80 \
  --alignment-mode fast --cascade-level fast --threads 8 --profile \
  --profile-json target/uniref50/profile.json

./target/release/kira-cluster build-index \
  target/uniref50/db target/uniref50/index \
  --k 14 --m 20 --pgm-epsilon 64

./target/release/kira-cluster search-index \
  target/uniref50/index tests/data/cluster_small.fasta target/uniref50/hits.tsv \
  --k 14 --m 20 --top-k 10 --min-seed-hits 2 \
  --prune-df-quantile 0.90 --max-seeds-per-query 128
```

For quicker comparisons, create a representative FASTA subset first and run:

```bash
KIRA_CLUSTER_BIN=target/release/kira-cluster \
  cargo bench --bench tool_compare -- data/uniref50_sample.fasta target/bench-uniref50-sample
```

## Sensitivity defaults

- Protein clustering seeds use a reduced (Dayhoff-6) amino-acid alphabet with a
  shorter default `k`, and alignment uses BLOSUM62 by default. This lets
  clustering surface homologs well below ~90% identity (exact long k-mers almost
  never collide across diverged proteins). `--sub-matrix simple` restores flat
  match/mismatch scoring.
- The default cascade is `fast` (positional, no gapped alignment) for
  predictable speed on any input size. Pass `--alignment-mode sensitive` (or
  `--cascade-level sensitive`) to enable banded, BLOSUM62-scored gapped alignment
  that also catches indel-containing homologs.
- `--cov-mode {0,1,2,3}` selects the MMseqs2-style coverage definition:
  `0` bidirectional (default), `1` target, `2` query, `3` target/query length ratio.

## Limitations / Experimental

- `serve` is a stub: it binds, answers a single request with a static JSON body,
  and exits. It is not a production server.
- `--seed-policy syncmer` is not implemented and emits `UNIMPLEMENTED_FLAG`.
- GPU (`--gpu --gpu-backend cuda`) is feature-gated (`--features cuda`) and falls
  back to CPU with a deterministic `KW2001 GPU_BACKEND_UNAVAILABLE_FALLBACK`
  warning when unavailable.
- `roaring` hot-posting storage is feature-gated (`--features roaring`).
- The index/search path uses exact k-mers (the reduced alphabet applies to
  clustering only).

## Notes

- If GPU is requested but unavailable, deterministic warning `KW2001 GPU_BACKEND_UNAVAILABLE_FALLBACK cuda` is emitted and CPU fallback is used.
- `benches/tool_compare.rs` runs optional MMseqs2/CD-HIT/VSEARCH comparisons when those tools are on `PATH`.
