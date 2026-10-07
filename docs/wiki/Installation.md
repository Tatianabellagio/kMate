# Installation

kMate is a Python package. It also calls external programs: **jellyfish** (k-mer
counting), **samtools** and **bcftools**.

## Install

```bash
mamba create -n kmate -c conda-forge -c bioconda kmate "python<3.13"
mamba activate kmate
```

`conda`/`micromamba` work the same way. This pulls kMate and the programs it calls,
`kmer-jellyfish`, `samtools` and `bcftools`.

The `"python<3.13"` pin matters only for kMate 0.1.4 and older, whose `build-index`
needed jellyfish's Python bindings (bioconda ships them for Python 3.9–3.12 only;
see [Troubleshooting](Troubleshooting)). From 0.1.5 kMate uses only the `jellyfish`
program and runs on any Python 3.9+; the pin is harmless.

### From source

To track the latest changes:

```bash
git clone https://github.com/Tatianabellagio/kMate.git
cd kMate
mamba create -n kmate -c conda-forge -c bioconda python numpy scipy pysam \
    kmer-jellyfish samtools bcftools
mamba activate kmate
pip install -e .
```

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

Use **0.1.5 or newer**. 0.1.0 cannot count k-mers; 0.1.2 cannot build a panel from a
Cactus pangenome without extra scripts, and its `--haploidize` empties haploid panels;
0.1.3 returns `NaN` without explanation for very-low-coverage pools; 0.1.3 and 0.1.4
need Python below 3.13 for `build-index`.

## From source

For development, or to track the latest changes:

```bash
git clone https://github.com/Tatianabellagio/kMate.git
cd kMate
mamba create -n kmate -c conda-forge -c bioconda python numpy scipy pysam \
    kmer-jellyfish samtools bcftools
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
