# kMate

Per-sample, per-record **allele-frequency estimation from pooled sequencing** against a
multi-founder reference panel. kMate runs a weighted k-mer Poisson EM on the 231-founder
simplex to estimate founder frequencies (`h`), then projects through a per-record presence/absence
matrix (`var_pa`, the founder × variant alt-allele matrix $V_\mathrm{pa}$) to allele frequencies for SNPs, indels, and SVs in a single pass.

The project root *is* the kMate estimator (renamed from the legacy `hapfire_sv/`).

## Start here

| Doc | What it is |
|---|---|
| [`HANDOFF.md`](HANDOFF.md) | Current session state + the production recipe (read first). |
| [`ALGORITHM.md`](ALGORITHM.md) | The kMate algorithm, math, and code wiring (code-verified source of truth). |
| [`BACKGROUND.md`](BACKGROUND.md) | Project framing, known biases, and design decisions. |
| [`SAVIO_HPC.md`](SAVIO_HPC.md) | Running on the Berkeley Savio cluster (partitions, sbatch recipes). |
| [`docs/`](docs/) | Analyses, investigations, methods writeups. `docs/PIPELINE_STATE.md` is the production-state source of truth. |

## Layout

```
src/          estimator code (em_solver, kmer_count, per_sample_per_chrom, kmer_pa builders); src/archive/ = retired variants
scripts/      helper/run scripts for the estimator; scripts/archive/ = superseded
tests/        dev/smoke tests + the production scale-out template (run_site_array_perchrom.sh)
data/         production kmer_pa matrices (kmer_pa_*, var_pa_*), sample lists, splits, small config

panel/        231-founder panel construction
  arch3/              arch decomposition → per-record var_pa matrices
  pangenie_index/     in-house k-mer index tables
  pangenie_genotyping/  short-read PanGenie genotyping of the 153 PG founders
  imputation/         deprecated (Beagle); kept for archaeology

benchmarks/   end-to-end accuracy benchmarks
  p80/   homogeneous 80-cactus-founder control
  p231/  full 231-founder headline benchmark
sims/         visor_freqk/ — shared pool-seq simulation framework

external/     vendored third-party tools (HapFIRE, pangenie-tools, genotyping-pipelines)
docs/         analyses + methods writeups
notebooks/    analysis notebooks
results/ plots/ logs/   outputs (gitignored)
archive/      retired exploration + one-off outputs
old_docs/     superseded docs (historical; not authoritative)
```

> `contamination_test/`, `preprocess_qc/`, `scratch/`, and `papers/` remain at the root as
> in-progress / working areas. The on-disk result-dir name `cactus_em_*` → `kmate_*` rename
> is a separate, deferred step (see `HANDOFF.md`).

## Production recipe (one-liner)

```bash
python src/per_sample_per_chrom.py \
    --kmer-pa-prefix data/kmer_pa_231_v3qc_v3_filt2inv/kmer_pa \
    --var-pa        panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz \
    --var-called panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz \
    --var-meta   panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz \
    --reads <r1.fq> <r2.fq> --sample <name> --out <out.tsv> \
    --threads 8 --chroms Chr1 --block-mode global --kmer-weight inv_mb
```

See `HANDOFF.md` for `global` vs `window` modes and the full production-state context.
