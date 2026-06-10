<p align="center">
  <img src="assets/kMate_logo.png" alt="kMate" width="420">
</p>

# kMate

Per-sample, per-record **allele-frequency estimation from pooled sequencing** against a
multi-founder reference panel. kMate runs a weighted k-mer Poisson EM on the 231-founder
simplex to estimate founder frequencies (`h`), then projects through a per-record presence/absence
matrix (`var_pa`, the founder × variant alt-allele matrix $V_\mathrm{pa}$) to allele frequencies for SNPs, indels, and SVs in a single pass.

The project root *is* the kMate estimator (renamed from the legacy `hapfire_sv/`).

## Install

kMate is pure Python — no build step. Create the `kmate` environment (mamba or conda) with its core dependencies, then run the estimator with that env's Python:

```bash
git clone https://github.com/Tatianabellagio/kMate.git
cd kMate
mamba create -n kmate -c conda-forge -c bioconda python numpy scipy pysam
```

Core deps: `numpy`, `scipy`, `pysam`. The authoritative cluster environment is documented in [`docs/PIPELINE_STATE.md`](docs/PIPELINE_STATE.md) §0.2; the run command is under **Production recipe** below.

## Start here

| Doc | What it is |
|---|---|
| [`docs/PIPELINE_STATE.md`](docs/PIPELINE_STATE.md) | **THE single source of truth — read first.** Production inputs (§0), run recipe (§0.1), environment (§0.2), what's deprecated. |
| [`ALGORITHM.md`](ALGORITHM.md) | The kMate algorithm, math, and code wiring (defers to PIPELINE_STATE for production paths). |
| [`BACKGROUND.md`](BACKGROUND.md) | Project framing, known biases, and design decisions. |
| [`SAVIO_HPC.md`](SAVIO_HPC.md) | Cluster ops (partitions, sbatch recipes). |
| [`HANDOFF.md`](HANDOFF.md) | Redirect stub — merged into `docs/PIPELINE_STATE.md` (2026-05-30). |

## Layout

The estimator itself lives in `src/` and is panel-agnostic; everything else is panel construction, the GrENE-Net scale-out, benchmarks, simulations, and supporting docs.

```
src/          estimator code (em_solver, kmer_count, per_sample_per_chrom, kmer_pa builders); src/archive/ = retired variants
scripts/      helper/run scripts for the estimator; scripts/archive/ = superseded
tests/        dev/smoke tests for the estimator (kmer-count, matrix builders, genome-wide validation); tests/archive/ = retired
grenenet/     GrENE-Net application of kMate (production scale-out runner for the evolved cohort); the estimator itself is panel-agnostic and lives in src/
data/         production kmer_pa matrices (kmer_pa_*, var_pa_*), sample lists, splits, small config

panel/        231-founder panel construction
  arch3/              arch decomposition → per-record var_pa matrices
  pangenie_index/     in-house k-mer index tables
  pangenie_genotyping/  short-read PanGenie genotyping of the 153 PG founders
  imputation/         deprecated (Beagle); kept for archaeology

benchmarks/   end-to-end accuracy benchmarks
  p80/   homogeneous 80-cactus-founder control
  p231/  full 231-founder headline benchmark
sims/         scripts/ — pool-seq simulation framework (mosaic builder + AF truth); see sims/README.md

external/     vendored third-party tools (gitignored; provenance + reconstruct steps in docs/EXTERNAL_DEPENDENCIES.md)
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

Everything derives from the arch3 panel VCF `merged_231_chr{N}_final.vcf.gz` — see
`docs/PIPELINE_STATE.md` §0. Env: `kmate`.

```bash
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python src/per_sample_per_chrom.py \
    --kmer-pa-prefix data/kmer_pa_231_arch3_filt2inv/kmer_pa \
    --var-pa     panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz \
    --var-called panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz \
    --var-meta   panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz \
    --reads <r1.fq> <r2.fq> --sample <name> --out <out.tsv> \
    --threads 8 --chroms Chr1 --block-mode global --kmer-weight inv_mb
```

See `docs/PIPELINE_STATE.md` §0.1 for `global` vs `window` modes and the full production-state context.

## Citation

kMate underlies the allele-frequency analyses of the GrENE-Net outdoor evolution experiment. If you use kMate, please cite:

> Xing Wu, Tatiana Bellagio, Yunru Peng, Lucas Czech, Meixi Lin, *et al.* (2026). Rapid adaptation and extinction across climates in synchronized outdoor evolution experiments of *Arabidopsis thaliana*. *Science* **391**, eadz0777. https://doi.org/10.1126/science.adz0777
>
> Wu, Bellagio, Peng, Czech & Lin contributed equally (co-first authors).

Preprint: bioRxiv 2025.05.28.654549 — https://doi.org/10.1101/2025.05.28.654549
