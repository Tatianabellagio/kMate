# Getting started — install kMate and run it on your own pool-seq

From zero to a per-record allele-frequency table. Two situations:

- **You already have a panel** (someone hands you `kmer_pa` + `var_pa` matrices) → §1, §2, §4.
- **You are building a panel** from your own founder assemblies → §1, §3, §4.

For the algorithm and the math see [`../ALGORITHM.md`](../ALGORITHM.md); for what is
production vs deprecated see [`PIPELINE_STATE.md`](PIPELINE_STATE.md) §0.

---

## 1. Install

```bash
mamba create -n kmate -c conda-forge -c bioconda kmate
mamba activate kmate
kmate selftest          # must print PASS before you go further
```

> ⚠️ **Check `kmate --version` reports 0.1.1 or newer.** The 0.1.0 build currently on
> bioconda is a pre-July-2026 snapshot with two problems: it depends on `jellyfish`
> instead of `kmer-jellyfish` (so it cannot count k-mers at all), and it predates the
> per-founder normalization fix — it has no `--normalize`, `--unit` or `--emit-af-se`,
> and would silently run the old EM in which k-mer-poor founders collapse toward zero.
> Until 0.1.1 is published, install from source (below).

`kmate selftest` runs a bundled tiny fixture end-to-end (k-mer counting → EM → AF
projection). It needs no network and no data of yours. **If it does not pass, stop** —
everything below depends on it.

> If you install the dependencies by hand, the k-mer counter is **`kmer-jellyfish`**.
> Plain `jellyfish` on conda-forge is an unrelated Python string-similarity library and
> leaves you with `jellyfish: command not found`.

## 2. What kMate needs

kMate estimates, for one pooled sample at a time, the **founder mixture** `h` and then
projects it to a frequency for every variant in the panel.

| input | what it is |
|---|---|
| pooled reads | paired FASTQ from one pool |
| `kmer_pa` | founder × k-mer presence/absence — the EM evidence |
| `var_pa` + `var_called` | founder × variant alt-allele, and the genotype-call mask |
| `meta` | per-record chrom/pos/ref_len/alt_len |

The matrices are **per chromosome** and are built once per panel (§3). They are large
(tens of GB for a 231-founder panel) and are not in the git repo.

## 3. Building a panel from founder assemblies

Only if you do not already have matrices. The full worked example — 8 *A. thaliana*
founders, every command, and the decisions that are panel-size dependent — is the MAGIC
runbook (see the project directory referenced in `PIPELINE_STATE.md` §6). In outline:

```
assemblies ──cactus-pangenome──▶ graph VCF ──build_kmers_tsv.py──▶ k-mer index
                                     │                                  │
                                     ├──── kmate build-var-pa ──▶ var_pa / var_called / meta
                                     └──── kmate build-kmer-pa ─▶ kmer_pa
```

```bash
# 1. pangenome (Minigraph-Cactus) -> a vcfbub-filtered VCF with real genotypes
cactus-pangenome <jobstore> <seqfile> --outDir out --outName mypanel \
    --reference <REFNAME> --haplo --vcf --gfa --gbz --giraffe

# 2. per-bubble k-mer index
python panel/pangenie_index/scripts/build_kmers_tsv.py \
    --vcf out/mypanel.vcf.gz --ref REF.fa --out index/ours -k 31 --haploid

# 3. founder x k-mer matrix, per chromosome
kmate build-kmer-pa --kmers index/ours_Chr1_kmers.tsv.gz --vcf out/mypanel.vcf.gz \
    --ref REF.iupacN.fa --chrom Chr1 --out kmer_pa/kmer_pa_Chr1 \
    --treat-missing-as-n --filter-production --min-ac 1 --invariant-margin 1

# 4. founder x variant matrices, per chromosome
kmate build-var-pa --vcf out/mypanel.vcf.gz --chrom Chr1 --out var_pa/var_pa_Chr1
```

**Three things that bite people:**

1. **`--min-ac` depends on panel size.** `--min-ac 2` (drop k-mers carried by one
   founder) is right for a large panel, where singletons are mostly noise. On a *small*
   panel those singletons are your most informative columns — on 8 founders `--min-ac 2`
   removed ~46% of k-mers. Use `--min-ac 1` for small panels.
2. **Two different reference FASTAs.** The graph build and the k-mer index take the
   plain reference; `build-kmer-pa` takes the **IUPAC→N** version, because Cactus
   normalises ambiguity codes inside the graph and the bubble flanks must match.
3. **Drop monomorphic records** (`AC=0 || AC=AN`) from the panel VCF before building
   `var_pa`, so the matrix is segregating-only.

## 4. Run kMate on a sample

```bash
# count k-mers once, query per chromosome (faster, byte-identical)
kmate build-kmer-db --reads R1.fq.gz R2.fq.gz --out pool.jf --threads 8

kmate run \
    --kmer-pa-prefix kmer_pa/kmer_pa \
    --var-pa     var_pa/var_pa_Chr1.var_pa.npz \
    --var-called var_pa/var_pa_Chr1.var_called.npz \
    --var-meta   var_pa/var_pa_Chr1.meta.npz \
    --reads R1.fq.gz R2.fq.gz --kmer-db pool.jf \
    --sample MYSAMPLE --out MYSAMPLE_Chr1.tsv \
    --threads 8 --chroms Chr1 \
    --unit chrom --kmer-weight uniform --normalize per_founder
```

### Choosing `--unit` — the one decision that matters most

**Pick it from your population's biology, not from a default.**

| your pools are… | use | why |
|---|---|---|
| selfing / inbred / a founder (F0) mix | `--unit chrom` + `--kmer-weight uniform` | ancestry is ~constant along a chromosome, so pooling all its k-mers gives the best-determined `h` |
| recombinant (a few generations of outcrossing) | `--unit bp --window-bp 10000` + `--kmer-weight inv_mb` | ancestry is a mosaic; it must be fit locally |

Getting this wrong is not subtle. On a selfing benchmark, fitting per LD-block instead of
per chromosome moved AF-MAE from 0.0033 to 0.0080 and outliers from 0.001% to 0.625%,
because low-diversity regions (centromeres) cannot identify founders from local k-mers
alone. Conversely, on a recombinant pool a chromosome-wide fit averages away the mosaic
you are trying to measure. Background:
[`EM_UNIT_CHOICE_AND_NONIDENTIFIABILITY.md`](EM_UNIT_CHOICE_AND_NONIDENTIFIABILITY.md).

`--normalize per_founder` should stay on: it divides each founder's EM update by its own
k-mer content, without which k-mer-poor founders collapse toward zero
([`FOUNDER_NORMALIZATION_FIX.md`](FOUNDER_NORMALIZATION_FIX.md)).

### Before you trust the numbers

QC your libraries **first** and drop the failures. The signature of a dead library is
normal trim survival, *low* duplicate loss, and flat low coverage across all chromosomes —
too few reads, not bad ones. Estimating on a 0.2× library produces output that looks like
data. Read preprocessing recipe:
[`PIPELINE_FASTQ_PREPROCESSING.md`](PIPELINE_FASTQ_PREPROCESSING.md).

## 5. Output

One row per panel record:

```
chrom  pos  ref_len  alt_len  alt_freq  info  n_called  se
```

- `alt_freq` — estimated alt-allele frequency in the pool
- `n_called` — founders with a genotype call at that record (panel support; independent of `h`)
- `se` — Wald SE using `n_called` as effective sample size
- `info` — `h`-weighted called mass

There is **no variant-type column**; derive it from the lengths:
SNP = `ref_len==1 and alt_len==1`; indel = `ref_len!=alt_len and max(...)<50`;
SV = `max(ref_len,alt_len)>=50`.

Run the same panel across samples and **row order is identical**, so you can join
row-wise. Do *not* join on `chrom:pos` — a position can carry more than one record, and a
position key will match the wrong alt allele.

## 6. Scaling to many samples

`grenenet/run_site_array_perchrom.sh` is a SLURM array runner: one task per sample, k-mer
DB built once and queried per chromosome, resumable. Point it at any panel with
`KMER_PA_PREFIX` / `VAR_PA_DIR` / `VAR_PA_TAG` — don't fork it.

The one lesson that generalises beyond kMate: the per-sample scratch files are the
bottleneck. Put the transient k-mer DB on **node-local / RAM** storage, never the shared
filesystem — at high concurrency that alone was up to a 48× difference in query time.
See [`../grenenet/README.md`](../grenenet/README.md).
