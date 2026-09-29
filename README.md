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

> ⚠️ **Do not `conda install kmate` yet.** bioconda currently serves **0.1.0**, which is
> broken two ways: it depends on `jellyfish` (a Python string-similarity library that
> ships no `jellyfish` binary, so kMate cannot count k-mers), and it was built from a
> pre-July-2026 snapshot lacking `--normalize`, `--unit` and `--emit-af-se` — it would
> silently run an older estimator. The fix is in review:
> [bioconda-recipes#69722](https://github.com/bioconda/bioconda-recipes/pull/69722). Until it merges, use one of the two routes below.

**Current release (0.1.1), from PyPI into a conda environment:**

```bash
mamba create -n kmate -c conda-forge -c bioconda python numpy scipy pysam kmer-jellyfish samtools
mamba activate kmate
pip install kmate               # bioconda's 0.1.0 is broken; see above
kmate selftest                  # must print PASS
```

**From source (for development, and to track `master`):**

```bash
git clone https://github.com/Tatianabellagio/kMate.git
cd kMate
mamba create -n kmate -c conda-forge -c bioconda python numpy scipy pysam kmer-jellyfish samtools
mamba activate kmate
pip install -e .
kmate selftest
```

The conda step installs the non-Python tools kMate shells out to: **`kmer-jellyfish`**
(the k-mer counter — *not* `jellyfish`) and `samtools`. The Python deps are `numpy`,
`scipy`, `pysam`.

### Verify the install

```bash
kmate selftest
```

This runs the bundled tiny fixture (a real Chr1 panel slice + a simulated 5-founder pool) end-to-end — exercising the full k-mer-count → EM → AF-projection path through `jellyfish`/`samtools` — and checks that the planted founder mixture is recovered. It takes a few seconds, needs no network, and prints `PASS` on a correct install. Run this **before** pointing kMate at your own data.

## Usage

kMate processes **one pooled sample at a time, per chromosome**.

**You need**
- **Pooled reads**: paired FASTQ (`R1.fq R2.fq`) of one pool/sample.
- **A reference panel** encoded as per-chromosome matrices: `kmer_pa` (k-mer × founder presence/absence), `var_pa` (founder × variant alt-allele), and record `meta`. Built once from your founders' VCF — see **[`docs/BUILDING_A_PANEL.md`](docs/BUILDING_A_PANEL.md)**. The 231-founder *Arabidopsis thaliana* panel used by GrENE-Net is available on request; the matrix files are large and are not stored in the Git repo.

**Run**

```bash
kmate run \
    --kmer-pa-prefix data/kmer_pa_231_arch3_filt2inv/kmer_pa \
    --var-pa     panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz \
    --var-called panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz \
    --var-meta   panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz \
    --reads R1.fq R2.fq --sample MYSAMPLE --out MYSAMPLE.tsv \
    --threads 8 --chroms Chr1 --kmer-weight uniform      # --unit chrom is the default
```

(`kmate run --help` lists every flag. Existing scripts that call `python src/per_sample_per_chrom.py ...` still work via thin shims that forward to the package.)

**Estimation unit** (`--unit`) — there is **one estimator**; the "mode" is just the unit it fits. Each unit is fit locally: *haploblock-collapse → EM → project*, with no anchor prior, no cross-window smoothing, and no fallback.
- `--unit chrom` (**default**): one founder mixture per chromosome. The default and the **selfing / inbred / F0** (GrENE-Net) production estimator — robust on uniform and sparse panels; also the only unit that supports `--h-only` and `--emit-af-se`.
- `--unit ld` (`--ld-r2 0.1`): r²-LD blocks derived from the panel's own `var_pa` (CompleteLDPartition). A per-block option — it **collapses in low-diversity blocks** (e.g. the centromere), so it is **wrong for selfing pools**; use it only for recombinant pools where fine per-block resolution helps.
- `--unit bp` (`--window-bp N`): fixed-bp windows, for **recombinant** pools.
- `--unit tsv` (`--blocks-tsv PATH`): explicit block partition.

`--block-mode global|window` are kept as **deprecated aliases** (`global`→`--unit chrom`, `window`→`--unit bp`). Recommended weighting is `--kmer-weight uniform` — the per-founder M-step normalization (kMate's default) removes the panel-completeness imbalance at its source, so `--kmer-weight inv_mb` is redundant (superseded 2026-07-06; see [`docs/FOUNDER_NORMALIZATION_FIX.md`](docs/FOUNDER_NORMALIZATION_FIX.md)).

**Haploblock collapse** (all units, on by default): before each EM, kMate computes the distinct k-mer haplotypes the panel actually resolves over the unit and fits those `K_b ≤ 231` haplotypes rather than assuming all 231 founders are separately identifiable, splitting each haplotype's frequency equally back to its members. `--haploblock-eps` sets the merge tolerance (default `0` = exact k-mer-identical, an exact no-op when all founders are distinct). This is the *block → haploblock → EM* design (like HARP/hapFIRE); it chiefly matters for finer units, where a small block may carry only a handful of haplotypes. See [`ALGORITHM.md`](ALGORITHM.md) §4.4.

The M-step normalization defaults to `--normalize per_founder`; pass `--normalize global` only to reproduce legacy (pre-2026-07-06) runs.

**Output**: a per-record TSV, one row per panel variant (SNP / indel / SV):

| chrom | pos | ref_len | alt_len | alt_freq | info | n_called | se |
|---|---|---|---|---|---|---|---|

`alt_freq` is the estimated alternate-allele frequency in the pool; `n_called` and `se` carry support/uncertainty. (`--var-called` adds a per-record called-mask; `--kmer-db` lets you count k-mers once and query per-chrom instead of re-scanning reads; `--hash-size` tunes the Jellyfish hash, e.g. lower it to `100M` on memory-capped jobs.)

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
