# grenenet — running kMate at cohort scale

kMate (the **method**) lives in `src/`. This folder is the **scale-out layer**:
how you run that method across thousands of pool-seq samples on a shared-filesystem
HPC cluster. Nothing here changes the math — it's all orchestration, I/O placement,
and SLURM plumbing. A single-sample / laptop user never needs any of it; they call
`src/per_sample_per_chrom.py` directly.

---

## The two layers (important mental model)

| Layer | Where | What changes it |
|---|---|---|
| **Method** — EM, AF projection | `src/per_sample_per_chrom.py` | `--block-mode`, `--kmer-weight`, the matrices |
| **Scale-out** — fan across N samples | `grenenet/*.sh` (this folder) | concurrency, memory, **DB I/O placement** |

The performance work below is all the **scale-out** layer. It is *deployment tuning*,
not part of kMate. Tell a new kMate user about the method knobs; only mention the
scale-out knobs if they're fanning thousands of samples on a shared filesystem.

---

## Single-phase runner (recommended): `run_site_array_perchrom.sh`

One SLURM array task per sample. Each task: builds a k-mer database once, then runs
the driver per chromosome (querying that DB), then concatenates to `OUT_DIR/<sample>.tsv`.
Per-chrom TSVs and the final TSV are skipped if present → preemption/requeue-safe.

```bash
sbatch --array=1-N%C --mem=32G \
  --export=ALL,MANIFEST=<tsv>,OUT_DIR=<dir>,BLOCK_MODE=global \
  grenenet/run_site_array_perchrom.sh
```
`MANIFEST` = TSV with header + `sample_id, reads_path[, reads_path2]`.

### Two performance levers (both automatic; here's *why* they exist)

**1. Count once (`--kmer-db`).** k-mer counting (scanning the read pool) dominated
runtime and was repeated once per chromosome (5× redundant). The runner now scans
the reads **once** into a database and *queries* it per chromosome. ~1.7–2× faster,
byte-identical output. This is genuinely useful even at small scale.

**2. Where the database lives (`JF_DIR`) — the shared-filesystem lever.** The
per-sample DB is **unique and cold** (2–3 GB, freshly written). When many tasks run
at once, all of them reading their own DBs off the *same shared filesystem* is the
dominant stall — measured **up to 48× query slowdown** at 75-way concurrency (worst
query 32 min vs ~40 s solo). Fix: put the DB on **fast node-local storage** so it
never touches the shared FS. The runner auto-picks:
`$SLURM_TMPDIR` (local SSD) → `/dev/shm` (RAM, if >4 GB free) → `OUT_DIR` (scratch, fallback).
A/B test (identical concurrency, only DB location changed): query-time **tail 1944 s → 39 s**, CPU efficiency ~27% → ~75%.
Override with `JF_DIR=<path>`. If you use `/dev/shm`, keep `--mem` ≥ ~28 GB (≈21 GB
peak compute + ~3 GB DB); the exit trap frees it even on preemption.

> **The transferable lesson** (true of *any* per-sample tool fanned across an HPC
> cluster, not just kMate): the per-sample scratch files are the bottleneck — put
> them on node-local / RAM storage, never the shared filesystem.

### Sizing (from clean, uncontended measurements)
- Peak RAM ≈ **21 GB** (global) → `--mem=32G` is safe; 64 GB is ~3× oversized.
- Per-sample wall: global ≈ 7–14 min, window ≈ 19 min (heavier block-EM).
- Cores: keep `--cpus-per-task=8` (window's block-EM uses up to 8 threads).

### Env knobs
`BLOCK_MODE` (**global**, default; `window` only if explicitly set) ·
`JF_DIR` (auto; override to force DB location) ·
`KMER_WEIGHT` (`inv_mb`) · `CHROMS` · `KMER_PA_PREFIX` · `VAR_PA_DIR`/`VAR_PA_TAG` ·
`OFFSET` (see chunking).

### Chunking past MaxArraySize (=1001 here)
No array can exceed 1000 indices, so a >1000-sample cohort is submitted in chunks
using `OFFSET` (manifest row = `OFFSET + task_id`):
```bash
for off in 0 1000 2000; do
  n=$(( TOTAL - off )); (( n > 1000 )) && n=1000; (( n <= 0 )) && break
  sbatch --array=1-${n}%C --mem=32G \
    --export=ALL,MANIFEST=<tsv>,OUT_DIR=<dir>,BLOCK_MODE=global,OFFSET=$off \
    grenenet/run_site_array_perchrom.sh
done
```

---

## Two-phase runner (optional): `run_phaseA/B/C` + `submit_two_phase_cohort.sh`

Phase A builds all DBs, Phase B is a flat `N×5` (sample,chrom) array, Phase C
concatenates. Maximizes fan-out for huge cohorts.

**Caveat that usually makes single-phase the better choice:** in two-phase the DB is
built by one job and read by 5 separate jobs that may land on *different* nodes — so
it **cannot** be node-local, putting Vessel-B I/O back on the shared filesystem. The
single-phase design keeps each sample's 5 chroms on one node, which is exactly what
lets the DB live in that node's RAM. Prefer single-phase unless you specifically need
the extra fan-out.

---

## Manifests (this cohort)
- `data/sample_manifest_usesample.tsv` — **2,168** samples = `usesample=True` = Table_S5
  = the hapFIRE analysis set ⭐ (recommended).
- `data/sample_manifest_savio.tsv` — 2,414 (all valid reads, incl. 246 QC-flagged).
- Reads are the **pooled** per-biological-sample dedup FASTQs (resequencing already
  merged). See the `grenenet-cohort-sample-structure` memory for the full decode.

See `docs/PIPELINE_STATE.md` §0/§0.1 for the production recipe.
