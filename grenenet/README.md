# grenenet — GrENE-Net application of kMate

kMate (the **algorithm**) lives in `src/`. This folder is the GrENE-Net
**project-specific application**: running kMate across the evolved GrENE-Net
pool-seq cohort (~2,415 `MLFH*` libraries) on the production 231-founder panel.
Keeping it out of `src/` preserves the separation between the reusable estimator
and this particular dataset's run.

## Contents
- `run_site_array_perchrom.sh` — SLURM array launcher (one task per sample) that
  runs `src/per_sample_per_chrom.py` per sample → per-record AF TSV. Defaults to
  the current production recipe: v3qc_v3 `kmer_pa`, arch3 `var_pa`/`var_called`,
  `--kmer-weight inv_mb`, window mode (10 kb). Reads a `MANIFEST` TSV
  (`sample_id, reads_path[, reads_path2]`); writes `OUT_DIR/<sample>.tsv`.

## Status (2026-05-30)
Updated to the current production recipe but **not yet run** — executing the full
cohort is a future session's task. Before running:
- **kmer_pa**: production filter is `filt2inv` (`ALGORITHM.md` §2.1); only `_filt2`
  is built on disk — rebuild and repoint `KMER_PA_PREFIX` first.
- **var_pa**: arch3 is Chr1 only — extend to Chr2-5 and set `CHROMS` accordingly.
- Evolved GrENE-Net reads are not PCR-free; ensure upstream dedup (clumpify) ran
  (`ALGORITHM.md` §10 M5, `docs/PIPELINE_FASTQ_PREPROCESSING.md`).

See `HANDOFF.md` for the production recipe and `docs/PIPELINE_STATE.md` §3.
