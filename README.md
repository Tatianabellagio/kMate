<p align="left">
  <img src="assets/kMate_logo.png" alt="kMate" width="420">
</p>

# kMate

[![license: MIT](https://img.shields.io/github/license/Tatianabellagio/kMate)](LICENSE)

Per-sample, per-record **allele-frequency estimation from pooled sequencing** against a
multi-founder reference panel. kMate runs a weighted k-mer Poisson EM on the founder
simplex to estimate founder frequencies (`h`), then projects through a per-record
presence/absence matrix (`var_pa`, the founder × variant alt-allele matrix $V_\mathrm{pa}$)
to allele frequencies for **SNPs, indels, and SVs in a single pass**, with no per-variant genotyping.

## Overview

<p align="center">
  <a href="poster_PEQG/poster_peqg.pdf">
    <img src="assets/poster_peqg.png" alt="kMate PEQG 2026 poster: tracking structural-variant trajectories across climates with alignment-free allele-frequency estimation" width="900">
  </a>
</p>

The picture above (our [PEQG 2026 poster](poster_PEQG/poster_peqg.pdf), click to enlarge) walks through the whole idea: the
GrENE-Net experiment evolved an equal mixture of **231 *Arabidopsis* founders** at 43 climate
sites over 3 years, pool-sequencing the surviving populations each generation. kMate takes those
pooled k-mer counts, solves a Poisson EM for the founder mixture against the panel's
`kmer_pa`/`var_pa` matrices, and reads out per-record allele frequencies for **SNPs *and* SVs** at
once. Benchmarked against simulated pools at 10× coverage, estimates track the truth closely, letting
us follow structural-variant frequency trajectories across climates (e.g. a 181-bp insertion in
the cold-regulated *COR413-PM2* gene, rising in cold gardens and falling in warm ones).

## Install

kMate has no build step. Create the `kmate` environment (mamba or conda) with its dependencies:

```bash
git clone https://github.com/Tatianabellagio/kMate.git
cd kMate
mamba create -n kmate -c conda-forge -c bioconda python numpy scipy pysam jellyfish samtools
```

Python deps: `numpy`, `scipy`, `pysam`. kMate also calls `jellyfish` (k-mer counting) and `samtools` (read handling), both installed by the command above.

## Usage

kMate processes **one pooled sample at a time, per chromosome**.

**You need**
- **Pooled reads**: paired FASTQ (`R1.fq R2.fq`) of one pool/sample.
- **A reference panel** encoded as per-chromosome matrices: `kmer_pa` (k-mer × founder presence/absence), `var_pa` (founder × variant alt-allele), and record `meta`. Built once from your founders' phased VCF (see [Building a panel](#building-a-panel)). The 231-founder *Arabidopsis thaliana* panel used by GrENE-Net is available on request; the matrix files are large and are not stored in the Git repo.

**Run**

```bash
python src/per_sample_per_chrom.py \
    --kmer-pa-prefix data/kmer_pa_231_arch3_filt2inv/kmer_pa \
    --var-pa     panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz \
    --var-called panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz \
    --var-meta   panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz \
    --reads R1.fq R2.fq --sample MYSAMPLE --out MYSAMPLE.tsv \
    --threads 8 --chroms Chr1 --kmer-weight inv_mb --block-mode global
```

**Estimator mode** (`--block-mode`)
- `global`: one founder mixture per chromosome. Use for **selfing / inbred / founder (F0)** pools.
- `window`: per-window mixture with HMM smoothing, for **recombinant** pools. `--block-mode window` alone reproduces the production "star2" recipe (10 kb windows, 5 smoothing passes).

**Output**: a per-record TSV, one row per panel variant (SNP / indel / SV):

| chrom | pos | ref_len | alt_len | alt_freq | info | n_called | se |
|---|---|---|---|---|---|---|---|

`alt_freq` is the estimated alternate-allele frequency in the pool; `n_called` and `se` carry support/uncertainty. (`--var-called` adds a per-record called-mask; `--kmer-db` lets you count k-mers once and query per-chrom instead of re-scanning reads.)

## How it works

From the read k-mer spectrum, kMate solves a weighted **Poisson EM on the 231-founder
simplex** for the founder mixture `h`, then projects `h` through the panel's
presence/absence matrix `var_pa` to a per-record allele frequency for SNPs, indels and SVs
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
grenenet/    GrENE-Net application: production scale-out over the evolved cohort
benchmarks/  end-to-end accuracy benchmarks (p80 control, p231 headline)
sims/        pool-seq simulation framework (AF truth); see sims/README.md
docs/        methods + analysis writeups
```

## Documentation

| Doc | What it is |
|---|---|
| [`docs/PIPELINE_STATE.md`](docs/PIPELINE_STATE.md) | Production inputs, run recipe, and environment; the project source of truth. |
| [`ALGORITHM.md`](ALGORITHM.md) | The kMate algorithm, math, and code wiring. |
| [`BACKGROUND.md`](BACKGROUND.md) | Project framing, known biases, and design decisions. |
| [`SAVIO_HPC.md`](SAVIO_HPC.md) | Cluster ops (partitions, sbatch recipes). |

## Using kMate

kMate is not yet published as a standalone method. If you are interested in using kMate
for your project, or in collaborating, please get in touch:

**Tatiana Bellagio** (tatianabellagio@gmail.com)

kMate was developed for, and underlies the allele-frequency analyses of, the GrENE-Net
outdoor evolution experiment in *Arabidopsis thaliana*.

## License

Released under the [MIT License](LICENSE).
