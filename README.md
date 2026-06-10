<p align="center">
  <img src="assets/kMate_logo.png" alt="kMate" width="420">
</p>

# kMate

Per-sample, per-record **allele-frequency estimation from pooled sequencing** against a
multi-founder reference panel. kMate runs a weighted k-mer Poisson EM on the founder
simplex to estimate founder frequencies (`h`), then projects through a per-record
presence/absence matrix (`var_pa`, the founder × variant alt-allele matrix $V_\mathrm{pa}$)
to allele frequencies for **SNPs, indels, and SVs in a single pass** — no per-variant genotyping.

## Install

kMate is pure Python — no build step. Create the `kmate` environment (mamba or conda) with its core dependencies:

```bash
git clone https://github.com/Tatianabellagio/kMate.git
cd kMate
mamba create -n kmate -c conda-forge -c bioconda python numpy scipy pysam
```

Core deps: `numpy`, `scipy`, `pysam`.

## Usage

kMate processes **one pooled sample at a time, per chromosome**.

**You need**
- **Pooled reads** — paired FASTQ (`R1.fq R2.fq`) of one pool/sample.
- **A reference panel** encoded as per-chromosome matrices: `kmer_pa` (k-mer × founder presence/absence), `var_pa` (founder × variant alt-allele), and record `meta`. Built once from your founders' phased VCF — see [Building a panel](#building-a-panel). This repo ships the 231-founder *Arabidopsis thaliana* panel.

**Run**

```bash
python src/per_sample_per_chrom.py \
    --kmer-pa-prefix panel/kmer_pa \
    --var-pa   panel/chr1/var_pa_chr1.var_pa.npz \
    --var-meta panel/chr1/var_pa_chr1.meta.npz \
    --reads R1.fq R2.fq --sample MYSAMPLE --out MYSAMPLE.tsv \
    --threads 8 --chroms Chr1 --block-mode global
```

**Estimator mode** (`--block-mode`)
- `global` — one founder mixture per chromosome. Use for **selfing / inbred / founder (F0)** pools.
- `window` — per-window mixture with HMM smoothing, for **recombinant** pools. `--block-mode window` alone reproduces the production "star2" recipe (10 kb windows, 5 smoothing passes).

**Output** — a per-record TSV, one row per panel variant (SNP / indel / SV):

| chrom | pos | ref_len | alt_len | alt_freq | info | n_called | se |
|---|---|---|---|---|---|---|---|

`alt_freq` is the estimated alternate-allele frequency in the pool; `n_called` and `se` carry support/uncertainty. (`--var-called` adds a per-record called-mask; `--kmer-db` lets you count k-mers once and query per-chrom instead of re-scanning reads.)

## How it works

From the read k-mer spectrum, kMate solves a weighted **Poisson EM on the 231-founder
simplex** for the founder mixture `h`, then projects `h` through the panel's
presence/absence matrix `var_pa` to a per-record allele frequency — SNPs, indels and SVs
together, in one pass. Processing one chromosome at a time keeps peak memory ~5× below a
genome-wide solve (per-chrom `h` agrees to ~0.1%). Full math + code wiring: [`ALGORITHM.md`](ALGORITHM.md).

## Building a panel

To run kMate on your own founder set you build the panel matrices once from a multi-founder
**phased VCF**: `var_pa` from the founder genotypes and `kmer_pa` from a k-mer index of the
founders. The builders live in [`panel/`](panel/) and [`data/`](data/). The bundled
231-founder *Arabidopsis* panel (used by GrENE-Net) and its exact construction are documented
in [`docs/PIPELINE_STATE.md`](docs/PIPELINE_STATE.md) §0.

## Repository layout

```
src/         the kMate estimator (em_solver, kmer_count, block_em, per_sample_per_chrom)
panel/       founder-panel construction (var_pa builders, k-mer index)
data/        prebuilt panel matrices (kmer_pa_*, var_pa_*) + sample lists
grenenet/    GrENE-Net application — production scale-out over the evolved cohort
benchmarks/  end-to-end accuracy benchmarks (p80 control, p231 headline)
sims/        pool-seq simulation framework (AF truth); see sims/README.md
docs/        methods + analysis writeups
```

## Documentation

| Doc | What it is |
|---|---|
| [`docs/PIPELINE_STATE.md`](docs/PIPELINE_STATE.md) | Production inputs, run recipe, and environment — the project source of truth. |
| [`ALGORITHM.md`](ALGORITHM.md) | The kMate algorithm, math, and code wiring. |
| [`BACKGROUND.md`](BACKGROUND.md) | Project framing, known biases, and design decisions. |
| [`SAVIO_HPC.md`](SAVIO_HPC.md) | Cluster ops (partitions, sbatch recipes). |

## Citation

kMate underlies the allele-frequency analyses of the GrENE-Net outdoor evolution experiment. If you use kMate, please cite:

> Xing Wu, Tatiana Bellagio, Yunru Peng, Lucas Czech, Meixi Lin, *et al.* (2026). Rapid adaptation and extinction across climates in synchronized outdoor evolution experiments of *Arabidopsis thaliana*. *Science* **391**, eadz0777. https://doi.org/10.1126/science.adz0777
>
> Wu, Bellagio, Peng, Czech & Lin contributed equally (co-first authors).

Preprint: bioRxiv 2025.05.28.654549 — https://doi.org/10.1101/2025.05.28.654549
