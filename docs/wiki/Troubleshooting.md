# Troubleshooting

Real errors kMate produces, and what they mean.

---

### `jellyfish: command not found`

The k-mer counter is missing. On conda-forge, `jellyfish` is an unrelated Python
library. The program you need is **`kmer-jellyfish`**:

```bash
mamba install -c conda-forge -c bioconda kmer-jellyfish
```

kMate 0.1.0 on bioconda depended on the wrong one; use 0.1.2 or newer (see
[Installation](Installation)).

---

### `WARNING: estimated coverage ... is below 1×` / `only N panel k-mers were observed`

The pool has too little sequence matching the panel for a stable estimate. Check the
library's depth first (a failed library looks like flat, very low coverage on every
chromosome), then that the reads come from a pool of this panel's founders. The
estimates are still written, but treat them as noise.

### `WARNING: no panel k-mer was observed in the reads`

Nothing in the reads matched the panel on that chromosome, so the mixture and every
allele frequency there are written as `NaN`. Usually the reads are not from a pool of
this panel's founders, or the library failed.

---

### `ModuleNotFoundError: No module named 'dna_jellyfish'`

`kmate build-index` uses jellyfish's Python bindings, and your environment's Python has
none. bioconda's `kmer-jellyfish` ships bindings only for Python 3.9–3.12; with a newer
Python (a fresh solve picks the newest) it installs an old build whose bindings are for
3.10 only. Recreate the environment with Python capped:

```bash
mamba create -n kmate -c conda-forge -c bioconda kmate "python<3.13"
```

Check with `python -c "import dna_jellyfish"`.

---

### `Panel VCF must be haploid (got len(GT)=2 ...)`

Your VCF has diploid genotypes. A kMate panel column is a *haplotype* (see
[Haplotypes and windows](Haplotypes-and-windows)), so one allele per genotype is required:

```bash
kmate decompose --haploidize --het missing ...   # inbred founders
kmate decompose --haploidize --het split   ...   # outbred, phased founders
```

See [Building a panel](Building-a-panel#some-founders-genotyped-from-short-reads).

---

### `Panel VCF must be biallelic (one ALT per record)`

Your VCF has records with several ALTs. Decompose with `kmate decompose`. Do **not** use
`bcftools norm -m -any`, which silently loses carriers on a pangenome graph.

---

### `--het split needs PHASED genotypes`

You asked to split samples into haplotypes, but a genotype is unphased (`0/1`). kMate
refuses rather than guess, because guessing invents haplotypes. Either phase the VCF,
or use `--het missing` if your founders are inbred lines.

---

### `kmate decompose --gfa` fails with `invalid literal for int()`

`annotate_vcf.py` needs **phased or haploid** genotypes (`0|1`, `1`), and stops on an
unphased one (`0/1`). A Minigraph-Cactus VCF is always phased, so this usually means the
VCF you passed is not the graph's own VCF. Use the `--annotated-vcf` /
`--biallelic-catalog` route for VCFs genotyped on the graph by another tool.

---

### `kmate selftest` fails

Something is wrong with the install itself; fix this before anything else. The usual
cause is the `jellyfish` mix-up above. Check `kmate --version` is 0.1.2 or newer.

---

### The job runs out of memory

The k-mer matrix is loaded densely: roughly **20 GB** for a few hundred founder haplotypes on
one chromosome. Request ~32 GB and run one chromosome at a time.

---

### Many samples at once are far slower than one

Almost always the k-mer databases sitting on shared network storage. Put them on
node-local disk or RAM (`$TMPDIR`, `/dev/shm`). We measured queries up to **48×**
slower when many jobs read their databases off a shared filesystem simultaneously.

---

### Everything runs, but the frequencies look wrong

Check these in order:

1. **Coverage.** Below ~1× the estimates are mostly noise. QC and drop failed
   libraries before estimating.
2. **`--unit`.** Selfing/inbred pools need `--unit chrom`; recombinant pools need
   `--unit bp`. The wrong one produces plausible-looking but wrong numbers; see
   [Running kMate](Running-kMate#--unit).
3. **The founder-haplotype mixture** (`*.h_per_chrom.npz`). If it collapsed onto a few
   haplotypes when you expected many, the problem is upstream of the frequencies.
4. **How you joined samples.** Merging on `chrom`/`pos` alone matches the wrong allele
   where a position carries more than one variant. Join row-wise.
5. **`n_called`.** Low-support variants are weak regardless of what `se` says.

---

### Still stuck

Open an issue at <https://github.com/Tatianabellagio/kMate/issues> with the command
you ran, the full error, and `kmate --version`.
