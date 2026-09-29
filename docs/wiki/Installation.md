# Installation

kMate is a Python package. It also calls two external programs: **jellyfish**
(k-mer counting) and **samtools**.

## Install

```bash
mamba create -n kmate -c conda-forge -c bioconda python numpy scipy pysam kmer-jellyfish samtools
mamba activate kmate
pip install kmate
```

`conda`/`micromamba` work the same way.

> ### Do not use `conda install kmate` yet
> The bioconda package is currently stuck at an old, broken 0.1.0: it pulls in the
> wrong `jellyfish` and cannot count k-mers at all. The fix is in review:
> [bioconda-recipes#69722](https://github.com/bioconda/bioconda-recipes/pull/69722). Install as shown above until it merges.

> ### The k-mer counter is `kmer-jellyfish`, not `jellyfish`
> On conda-forge, `jellyfish` is an unrelated Python string-similarity library that
> installs no `jellyfish` program. If you install that by mistake you will get
> `jellyfish: command not found` when you run kMate.

## Check it works

```bash
kmate selftest
```

This runs a tiny bundled example end to end: k-mer counting, the EM, and the
frequency projection.

It should end with:

```
PASS: kMate is correctly installed and working.
```

**If it does not print PASS, stop and fix that first.** Everything else depends on it.

## Check the version

```bash
kmate --version
```

Use **0.1.2 or newer**. Older builds lack `--normalize`, `--unit` and `--emit-af-se`
and silently use a superseded estimator.

## From source

For development, or to track the latest changes:

```bash
git clone https://github.com/Tatianabellagio/kMate.git
cd kMate
mamba create -n kmate -c conda-forge -c bioconda python numpy scipy pysam kmer-jellyfish samtools
mamba activate kmate
pip install -e .
kmate selftest
```

## What you get

```
kmate run             estimate frequencies for one pool
kmate build-index     build the k-mer index for a panel
kmate build-kmer-pa   build the founder x k-mer matrix
kmate build-var-pa    build the founder x variant matrix
kmate decompose       make a multi-allelic VCF biallelic
kmate build-kmer-db   count a pool's k-mers once, reuse per chromosome
kmate selftest        verify the install
```

Run `kmate <command> --help` for the options of any one of them.
