<p align="left">
  <img src="assets/kMate_logo.png" alt="kMate" width="420">
</p>

# kMate

[![license: MIT](https://img.shields.io/github/license/Tatianabellagio/kMate)](LICENSE)

kMate estimates **allele frequencies for SNPs, short indels, and large structural variants**
from pooled sequencing of populations descended from a known set of founders. It follows the
logic of HARP ([Kessner et al. 2013](https://doi.org/10.1093/molbev/mst016)) but replaces
HARP's per-base read likelihood with a k-mer–count likelihood, making it **alignment-free and
pangenome-native**. kMate works in two stages: it first infers the founder mixture `h` from
observed k-mer counts by Poisson EM, then projects `ĥ` through a precomputed founder × variant
matrix to obtain an allele frequency for every record **in a single pass**.

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

```bash
mamba create -n kmate -c conda-forge -c bioconda kmate
mamba activate kmate
kmate selftest
```

This installs kMate and the programs it calls, **`kmer-jellyfish`** (the k-mer counter)
and `samtools`.

### From source

```bash
git clone https://github.com/Tatianabellagio/kMate.git
cd kMate
mamba create -n kmate -c conda-forge -c bioconda python numpy scipy pysam kmer-jellyfish samtools
mamba activate kmate
pip install -e .
kmate selftest
```

### Verify the install

```bash
kmate selftest
```

This runs the bundled tiny fixture (a real Chr1 panel slice + a simulated 5-founder pool) end-to-end — exercising the full k-mer-count → EM → AF-projection path through `jellyfish`/`samtools` — and checks that the planted founder mixture is recovered. It takes a few seconds and prints `PASS` on a correct install. Run this **before** pointing kMate at your own data.

## Usage

Two steps: turn your founder haplotypes into matrices once, then run each pool against
them.

### 1. Build the panel

One VCF of your founder haplotypes becomes two matrices, reused for every pool:

```bash
kmate build-index   --vcf panel.vcf.gz --ref REF.fa --out index/ours -k 31 --haploid

kmate build-kmer-pa --kmers index/ours_Chr1_kmers.tsv.gz \
                    --vcf panel.vcf.gz --ref REF.fa --chrom Chr1 \
                    --out kmer_pa/kmer_pa_Chr1 \
                    --treat-missing-as-n --filter-production

kmate build-var-pa  --vcf panel.vcf.gz --chrom Chr1 --out var_pa/var_pa_Chr1
```

`kmer_pa` is the founder-haplotype × k-mer matrix, the evidence the reads are compared
against; `var_pa` is the founder-haplotype × variant matrix the mixture is projected
through, written alongside `var_called` and `meta`.

The panel VCF must be **haploid, biallelic and sequence-resolved**.
**[`docs/BUILDING_A_PANEL.md`](docs/BUILDING_A_PANEL.md)** covers how to get there,
including building one from assemblies with `cactus-pangenome`, and `kmate decompose`
for multi-allelic input.

The 231-founder *Arabidopsis thaliana* panel used by GrENE-Net is available on request.

### 2. Run a pool

```bash
kmate run \
    --kmer-pa-prefix kmer_pa/kmer_pa \
    --var-pa-prefix  var_pa/var_pa \
    --reads R1.fq.gz R2.fq.gz \
    --sample MYSAMPLE --out MYSAMPLE_Chr1.tsv \
    --threads 8 --chroms Chr1
```

Both `--*-prefix` options are prefixes: kMate appends `_<CHROM>.kmer_pa.npz` and
`_<CHROM>.{var_pa,var_called,meta}.npz`. `kmate build-kmer-db` counts a pool's k-mers
once so `--kmer-db` can query them per chromosome.

### Choosing `--unit`

`--unit` sets the span over which one founder-haplotype mixture is estimated. Pick it
from your population's biology:

| your pools are | use | why |
|---|---|---|
| selfing, inbred, or a founder (F0) mix | `--unit chrom` (default) | ancestry is constant along a chromosome, so pooling all its k-mers determines the mixture best |
| recombinant | `--unit bp --window-bp 10000` | ancestry is a mosaic and must be fitted locally |

Fitting per LD block on a selfing pool moved AF error from 0.0033 to 0.0080 and outliers
from 0.001% to 0.625%, because low-diversity regions cannot distinguish founder
haplotypes from local k-mers alone. Windows are fitted independently by default;
`--smooth-windows` enables the anchor and cross-window smoothing.

`kmate run --help` lists every flag, and the
[Command reference](https://github.com/Tatianabellagio/kMate/wiki/Command-reference)
gives every default.

**Output**: a per-record TSV, one row per panel variant (SNP / indel / SV):

| chrom | pos | ref_len | alt_len | alt_freq | info | n_called | se |
|---|---|---|---|---|---|---|---|

`alt_freq` is the estimated alternate-allele frequency in the pool; `n_called` is how many
founder haplotypes had a genotype call there, and `se` is the binomial SE from it.

## How it works

From the read k-mer spectrum, kMate solves a weighted **Poisson EM on the 231-founder
simplex** for the founder mixture `h`, then projects `h` through the panel's
presence/absence matrix `var_pa` to a per-record allele frequency for SNPs, indels and SVs
together, in one pass. Processing one chromosome at a time keeps peak memory ~5× below a
genome-wide solve (per-chrom `h` agrees to ~0.1%). Full math + code wiring: [`ALGORITHM.md`](ALGORITHM.md).

## Building a panel

**→ [`docs/BUILDING_A_PANEL.md`](docs/BUILDING_A_PANEL.md)** — the step-by-step guide: what
the panel VCF must satisfy, how to produce one (from assemblies via a pangenome graph, from
an existing graph VCF, or from a phased callset), the three build commands, how to pick
`--min-ac` for your panel size, and how to check the result before you use it.

In short: you build the matrices once from a multi-founder **haploid, biallelic** VCF —
`var_pa` from the founder genotypes, `kmer_pa` from a k-mer index of the same VCF. The
231-founder *Arabidopsis* panel used by GrENE-Net is available on request (the matrices are
large and not in the repo); its exact construction is in
[`docs/PIPELINE_STATE.md`](docs/PIPELINE_STATE.md) §0.

## Repository layout

```
src/kmate/   the kMate package (em_solver, kmer_count, block_em, per_sample_per_chrom, cli, selftest)
pyproject.toml, conda/   packaging: pip-installable `kmate` CLI + conda recipe
panel/       founder-panel construction (var_pa builders, k-mer index)
data/        prebuilt panel matrices (kmer_pa_*, var_pa_*) + sample lists
grenenet/    GrENE-Net application: production scale-out over the evolved cohort
benchmarks/  end-to-end accuracy benchmarks (p80 control, p231 headline)
sims/        pool-seq simulation framework (AF truth); see sims/README.md
docs/        methods + analysis writeups
```


> **`archive/` trees are not tracked.** Retired work (`archive/`, `archive_gea/`, and
> `archive/` directories nested anywhere, e.g. `src/archive/`) stays on disk in a working
> checkout but is deliberately excluded from the repository. Links to `archive/…` paths in
> the docs therefore resolve only in a local checkout, not on GitHub.

## Documentation

| Doc | What it is |
|---|---|
| [`docs/GETTING_STARTED.md`](docs/GETTING_STARTED.md) | **Start here** — install, inputs, building a panel, running a sample, choosing `--unit`. |
| [`docs/PIPELINE_STATE.md`](docs/PIPELINE_STATE.md) | Production inputs, run recipe, and environment; the project source of truth. |
| [`ALGORITHM.md`](ALGORITHM.md) | The kMate algorithm, math, and code wiring. |
| [`BACKGROUND.md`](BACKGROUND.md) | Project framing, known biases, and design decisions. |
| [`SAVIO_HPC.md`](SAVIO_HPC.md) | Cluster ops (partitions, sbatch recipes). |

## Using kMate

kMate is not yet published as a standalone method. If you are interested in using kMate
for your project, or in collaborating, please get in touch:

**Tatiana Bellagio** (tatianabellagio@gmail.com)

kMate was developed for, and underlies the allele-frequency analyses of, the [GrENE-Net](https://www.science.org/doi/10.1126/science.adz0777)
outdoor evolution experiment in *Arabidopsis thaliana*.

## License

Released under the [MIT License](LICENSE).
