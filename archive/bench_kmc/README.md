# Archived: KMC 3 / rapidgzip k-mer counting benchmark

Bench scripts + SLURM logs from a 2026-05-02 test of replacing jellyfish in
the per-sample pipeline. **Verdict: jellyfish stays.**

**Canonical reasoning lives at the project level — read these instead:**
- `METHODS_TRIED.md` § 5 — "HARP REJECT" entry; pointer to next real lever
- `HANDOFF.md` 2026-05-02 entry — full headline table, CLI-bug postmortem,
  literature pointer (SSHash / LP-MPHF) for the next time someone wants to
  speed up k-mer counting

This directory is a candidate for deletion once disk is needed. Nothing here
is referenced by production code paths.

## What's in this directory

```
run_bench.sh              KMC 3 vs jellyfish driver
run_rapidgzip_bench.sh    rapidgzip v1 (had the CLI bug)
run_rapidgzip_bench2.sh   rapidgzip v2 (fixed CLI: cat R1 R2 | rapidgzip)
build_query_fa.py         Chr1 query FASTA builder
compare_outputs.py        byte-level / count-equality check across tools
bench_59298.{out,err}     KMC 3 SLURM logs (wall, RSS, kmer-equality)
rgz_bench_59320.{out,err} rapidgzip v1 SLURM logs (with the bug)
rgz_bench2_59551.{out,err} rapidgzip v2 SLURM logs (corrected)
```

The three SLURM workdirs (`work_59298` 8.5 GB, `rgz_59320` 17 GB,
`rgz2_59551` 11 GB) and the regeneratable `query_chr1.fa` (310 MB) were
deleted on 2026-05-08 to recover ~36 GB. They were intermediate read counts
+ k-mer databases, fully reproducible from the scripts above.

## How to re-run (if needed)

```bash
sbatch archive/bench_kmc/run_bench.sh             # KMC 3 vs jellyfish
sbatch archive/bench_kmc/run_rapidgzip_bench2.sh  # corrected rapidgzip
```
