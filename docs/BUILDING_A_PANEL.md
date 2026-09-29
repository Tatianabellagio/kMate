# Building a panel — from a VCF or pangenome to the kMate matrices

`kmate run` consumes four per-chromosome files. This page is how you produce them.

---

## The panel files

Built once per founder set, then reused for every sample.

| file | shape | what it is | built by |
|---|---|---|---|
| `kmer_pa_<CHR>.kmer_pa.npz` | founder × k-mer | which founders carry each panel k-mer — **the EM evidence** | `kmate build-kmer-pa` |
| `kmer_pa_<CHR>.meta.npz` | — | k-mer → bubble/chromosome map | same |
| `var_pa_<CHR>.var_pa.npz` | founder × variant | which founders carry the ALT — **the projection target** | `kmate build-var-pa` |
| `var_pa_<CHR>.var_called.npz` | founder × variant | which founders have a genotype call (not `./.`) | same |
| `var_pa_<CHR>.meta.npz` | — | per-record `chrom`, `pos`, `ref`, `alt` | same |

Both matrices are derived from **the same panel VCF**. Build them from different
VCFs and the EM and the projection disagree silently.

Producing them is three steps:

```
panel VCF ──(1) kmate build-index ────▶ k-mer index ──┐
    │                                                ├──(2) kmate build-kmer-pa ──▶ kmer_pa
    └────────────────────────────────────────────────┘
    └──(3) kmate build-var-pa ──▶ var_pa + var_called + meta
```


---

## Input requirements

### Input variants (the panel VCF)

- **Multi-sample.** One column per founder. Sample names become the founder axis, and
  that order is shared by both matrices — don't reorder between builds.
- **Haploid.** Exactly one allele per `GT` (`0`, `1`, `.` — not `0/1` or `0|1`). Both
  builders **abort with an error** on diploid genotypes rather than guess, because
  `kmer_pa` reconstructs one sequence per founder while `var_pa` uses a carrier rule;
  on a diploid GT those two would disagree. Haploidize before building.
- **Biallelic — one ALT per record.** ⚠️ **This is not checked, and violating it is
  silently wrong.** `var_pa` records `alts[0]` as *the* ALT but marks a founder as a
  carrier if any GT allele is `> 0` — so at a multi-allelic record a founder carrying
  ALT 2 is recorded as carrying ALT 1. `kmer_pa` handles multi-allelic records
  correctly, so the two matrices quietly stop agreeing. Decompose first.
- **Sequence-resolved.** Explicit REF/ALT sequences. Symbolic alleles (`<DEL>`, `<INS>`,
  `<CNV>`) cannot be turned into k-mers.
- **Non-overlapping.** One variant per locus per haplotype — for graph VCFs this means
  top-level bubbles only, no nested snarls (that is what `vcfbub` gives you).
- **bgzip-compressed and indexed** (`.vcf.gz` + `.tbi`/`.csi`). Both builders fetch by
  region. (The step-1 index builder also accepts a plain `.vcf`; its `--vcf` help text
  saying "uncompressed" is stale — gzip is auto-detected.)
- **Segregating-only, recommended.** Drop records where the ALT is carried by nobody or
  by every called founder (`AC=0 || AC=AN`). They cost memory and contribute nothing:
  an invariant k-mer adds the same constant to every founder's EM term.

```bash
bcftools view -e 'INFO/AC=0 || INFO/AC=INFO/AN' -Oz -o panel.seg.vcf.gz panel.vcf.gz
bcftools index -t panel.seg.vcf.gz
```

### Input reference

- **Indexed FASTA** (`.fai`) matching the VCF's coordinates and contig names.
- **One reference, used for every step.** IUPAC ambiguity codes are fine: kMate
  normalises any non-ACGTN base to `N` when it reads the reference, which is what
  Minigraph-Cactus does when building the graph, so the reference and the graph's REF
  alleles agree. k-mers spanning such a position contain `N` and are skipped.
  `build-kmer-pa` reports how many bases it normalised.

  *(Earlier versions of this pipeline required a separately prepared IUPAC→N copy of the
  reference for `build-kmer-pa` while other steps used the plain one — mixing them up
  silently corrupted bubble flanks. That is handled internally now; a single reference is
  correct everywhere. On TAIR10 the two files differed at exactly 469 positions, all
  IUPAC→N.)*

---

## How to generate the panel VCF

Three starting points.

### From assemblies, via a pangenome graph (recommended)

Gives the full SNP+indel+SV spectrum with real genotypes.

```bash
# seqfile: "<name><TAB><path>", reference FIRST
cactus-pangenome <jobstore> panel.seqfile \
    --outDir out --outName mypanel \
    --reference REFNAME --haplo --vcf --gfa --gbz --giraffe \
    --maxDisk 500G
```

`out/mypanel.vcf.gz` is already `vcfbub`-filtered (top-level bubbles). Genotypes are
haploid with explicit `0` / `1` / `.`, so it satisfies the requirements above — except
that it may be multi-allelic, so decompose (below).

Put the jobstore on **node-local disk**, not shared scratch. Pin **Cactus ≥ 3.1.0** if
your reference contains IUPAC codes: 3.1.0 normalises them to `N`, 2.9.3 errors out
(`Non-ACGTN character 'y'`).

### From an existing graph VCF

Reduce an existing graph VCF (HPRC, a published pangenome) to top-level bubbles:

```bash
vcfbub -l 0 -a 100000 --input graph.vcf.gz | bgzip > panel.vcf.gz
bcftools index -t panel.vcf.gz
```

### From a phased multi-sample callset (no graph)

Workable, but you only get what the callset contains — typically SNPs and short indels,
no SVs. You still need it biallelic and haploid.

## Making a multi-allelic VCF biallelic — `kmate decompose`

Graph VCFs are usually multi-allelic, and kMate requires one ALT per record.

**Do not decompose by realignment.** `bcftools norm -m -any` (and `--atomize`,
`vcfwave`) split each record by pairwise REF↔ALT alignment. On a pangenome graph that
scatters carriers across shifted positions and silently drops them where ALT paths
converge on the same atomic variant — measured here at up to **~99% carrier loss** at a
SNP co-located with a multi-allelic indel. It is silent: nothing errors.

The correct method is **symbolic-ID propagation**, developed for the Human Pangenome
Reference Consortium: each atomic variant nested in a bubble carries a symbolic ID, and
genotypes move from the multi-allelic record to the biallelic catalog by *matching IDs*,
never by alignment.

```bash
kmate decompose \
    --annotated-vcf      annotated_multiallelic.vcf.gz \   # from annotate_vcf.py
    --biallelic-catalog  annotated_biallelic.vcf.gz \      # from annotate_vcf.py
    --genotyped-vcf      your_genotyped.vcf.gz \
    --convert-to-biallelic /path/to/convert-to-biallelic.py \
    --haploidize \
    --out panel.vcf.gz
```

`--haploidize` also makes genotypes haploid, which the builders require — **a kMate panel
column is a haplotype**, because `build_kmer_pa` reconstructs exactly one sequence per
founder. That is a property of the model, not of any particular species. How to get there
depends on what your founders are:

| your founders are | use | what happens |
|---|---|---|
| **inbred lines** (*Arabidopsis* accessions, MAGIC/RIL founders, NAM parents) | `--het missing` (default) | `0/0`→`0`, `1/1`→`1`, **het→`.`**, `./.`→`.` |
| **outbred and phased** (HPRC-style assemblies, phased diploids) | `--het split` | each sample becomes two haplotype columns `sample.h1`/`sample.h2`; the founder axis doubles |

With inbred founders a heterozygous call is more likely a genotyping artefact than real
diploidy, so marking it missing lets the AF projection *exclude* that founder at that
record rather than invent a REF or ALT call (precedent for inbred *Arabidopsis* panels:
Arouisse et al. 2020, Plant J 102:872–882). The command reports the het→missing rate and
**warns if it exceeds 10%**, which usually means your founders are not what this policy
assumes.

With outbred founders that rule would discard half your data and bias AF, so `--het split`
keeps both haplotypes and kMate estimates a frequency per haplotype — which is what a
founder-mixture model means for an outbred panel. It requires **phased** genotypes; an
unphased heterozygote is an error rather than a coin flip, because guessing phase would
fabricate haplotypes that were never observed.

> ⚠️ `--het split` has so far only been exercised on synthetic input — see
> [`OPEN_ITEMS.md`](OPEN_ITEMS.md) §1 before relying on it for a real phased panel.
> In particular it assumes diploid, fully-phased genotypes; partially phased blocks,
> `PS` phase-set tags and ploidy ≠ 2 are not handled.

### Acknowledgement

The decomposition method is not kMate's — it is **HPRC symbolic-ID propagation**, run via
`annotate_vcf.py` (HPRC `prepare-vcf-MC`) and `convert-to-biallelic.py`
([eblerjana/pangenie-tools](https://github.com/eblerjana/pangenie-tools)). Neither is
bundled; obtain them separately. Only the `INFO/ID` transfer (`kmate transfer-id`) is
kMate's own — it exists because `bcftools annotate -c INFO/ID` corrupts the
angle-bracketed graph-node IDs.

If the decomposition matters to your results, cite Ebler et al. (2022) *Nat Genet*
54:518–525 and Liao et al. (2023) *Nature* 617:312–324. `kmate decompose --citation`
prints this.

**Inherited limitation:** `annotate_vcf.py` builds its atomic catalog with `vcfwave`
internally, so a small fraction of atomic variants (~0.8% on this project's panel) are
absent from the catalog and cannot be recovered by ID matching.

---

## Building the matrices

Per chromosome — loop or use a job array. `$CHR` is the VCF's contig name.

### Step 1 — k-mer index

```bash
kmate build-index \
    --vcf  panel.vcf.gz \
    --ref  REF.fa \
    --out  index/ours \
    -k 31 --haploid \
    --jellyfish-threads 8 --jellyfish-hash 3000000000
```

Writes `index/ours_<CHR>_kmers.tsv.gz`: one row per bubble with its unique k-mers.
Reproduces `PanGenie-index`'s output; full algorithm in
[`../panel/pangenie_index/kmer_index_construction.md`](../panel/pangenie_index/kmer_index_construction.md).
`--haploid` is required for haploid GTs.

### Step 2 — `kmer_pa`

```bash
kmate build-kmer-pa \
    --kmers index/ours_${CHR}_kmers.tsv.gz \
    --vcf   panel.vcf.gz \
    --ref   REF.fa \
    --chrom $CHR \
    --out   kmer_pa/kmer_pa_${CHR} \
    --treat-missing-as-n \
    --filter-production --min-ac 1 --invariant-margin 1
```

- `--treat-missing-as-n` — a `./.` genotype masks that founder's sequence over the
  record with `N`, dropping the overlapping k-mers. The alternative (treating `./.` as
  REF) invents a confident genotype and biases the panel toward reference.
- `--filter-production` — drop uninformative k-mer columns at build time (see below).

### Step 3 — `var_pa`

```bash
kmate build-var-pa --vcf panel.vcf.gz --chrom $CHR --out var_pa/var_pa_${CHR}
```

Writes `.var_pa.npz`, `.var_called.npz`, `.meta.npz`. The call mask is not optional
bookkeeping: `AF = (h·var_pa)/(h·var_called)` is what keeps `./.` from being counted as
REF and under-estimating AF wherever missingness is high.

---

## Choosing `--min-ac` — it depends on panel size

`--filter-production` keeps k-mer columns with `min_ac ≤ a_k ≤ F − invariant_margin`,
where `a_k` is how many founders carry that k-mer.

| founders | use | why |
|---|---|---|
| large (order 100+) | `--min-ac 2` | a k-mer in exactly one founder is mostly private repeat or sequencing noise, and gives no cross-founder discrimination |
| small (order ≲ 20) | `--min-ac 1` | founder-private k-mers are your *most* informative columns |

Measured on an 8-founder panel: `--min-ac 2` removed **~46% of k-mers, nearly all
`a_k == 1`**. The default is 2 — override it for a small panel.

`--invariant-margin 1` drops `a_k == F` (carried by everyone): `μ_k = 1` for every `h`,
so it contributes an identical constant to every founder and only dilutes the update.

---

## Check the panel before you use it

```python
import numpy as np, scipy.sparse as sp
K = sp.load_npz("kmer_pa/kmer_pa_Chr1.kmer_pa.npz")
V = sp.load_npz("var_pa/var_pa_Chr1.var_pa.npz")
C = sp.load_npz("var_pa/var_pa_Chr1.var_called.npz")

assert K.shape[0] == V.shape[0] == C.shape[0], "founder axis differs between matrices"
assert V.shape == C.shape

ac = np.asarray(K.sum(0)).ravel()          # founders per k-mer
print("founders:", K.shape[0], " k-mers:", K.shape[1], " records:", V.shape[1])
print("a_k range:", ac.min(), "-", ac.max(), "(expect >=min_ac and <=F-1)")
print("k-mers per founder: min", np.asarray(K.sum(1)).ravel().min(),
      " max", np.asarray(K.sum(1)).ravel().max())
print("records with zero called founders:",
      int((np.asarray(C.sum(0)).ravel() == 0).sum()), "(should be 0)")
```

What to look at:

- **Founder counts must match across all three matrices**, in the same order.
- **`a_k` range** must respect your filter — if `a_k` reaches `F` you did not apply
  `--invariant-margin`, if it hits 1 with `--min-ac 2` the filter did not run.
- **k-mers per founder**: a founder with far fewer than the others contributes less
  evidence. Large spreads are real (long-read assemblies realise more private k-mers
  than short-read-genotyped founders) and are exactly what `--normalize per_founder`
  corrects for — but a founder near zero usually means a broken sample column.
- **Records with zero called founders** make AF undefined; there should be none if you
  filtered to segregating-only.

Then run `kmate selftest` (bundled fixture, not your panel) to confirm the install
itself is sound, and estimate one known sample before launching a cohort.

---

## Layout the scale-out runner expects

`grenenet/run_site_array_perchrom.sh` builds paths as
`$VAR_PA_DIR/<chr-lowercase>/${VAR_PA_TAG}_<chr-lowercase>.{var_pa,var_called,meta}.npz`
and `${KMER_PA_PREFIX}_<CHR>.kmer_pa.npz`. So:

```
panel/
  kmer_pa/kmer_pa_Chr1.kmer_pa.npz        # KMER_PA_PREFIX=panel/kmer_pa/kmer_pa
  var_pa/chr1/var_pa_mypanel_chr1.var_pa.npz    # VAR_PA_DIR=panel/var_pa
  var_pa/chr1/var_pa_mypanel_chr1.var_called.npz  #  VAR_PA_TAG=var_pa_mypanel
  var_pa/chr1/var_pa_mypanel_chr1.meta.npz
```

Note the case difference — `kmer_pa` uses `Chr1`, the `var_pa` directories use `chr1`.
Symlinks are fine if you built them under one spelling.
