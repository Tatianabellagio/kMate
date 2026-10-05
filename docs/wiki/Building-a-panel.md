# Building a panel

Do this **once** per founder-haplotype set. Afterwards you reuse the matrices for every
pool.

```
founder assemblies ──▶ graph VCF ──────────────────▶ k-mer index ──┐
                           └──decompose──▶ panel VCF ──────────────┴──▶ kmer_pa
                                               └──────────────────────▶ var_pa + var_called + meta
```

Step 1 is `cactus-pangenome`; steps 2–4 are `kmate` commands, run **per chromosome**.

---

## Step 1: build a pangenome from your founder assemblies

The usual starting point is **long-read assemblies of your founders**. Build a
pangenome graph from them; this gives the full variant spectrum (SNPs, indels **and**
structural variants) with real genotypes.

```bash
cactus-pangenome <jobstore> panel.seqfile \
    --outDir out --outName mypanel \
    --reference REFNAME --haplo --vcf --gfa --gbz --giraffe
```

`panel.seqfile` is `<name><TAB><path>`, with the **reference first**, then one line
per **haplotype assembly**:

```
REFNAME     /path/to/reference.fa
founderA    /path/to/founderA.fa
founderB    /path/to/founderB.fa
```

- An **inbred** founder has one assembly: one line, one panel column.
- An **outbred** founder with a **haplotype-resolved** (phased) assembly has two: give
  each its own name (`plantA_hap1`, `plantA_hap2`) and you get two panel columns.

Either way `out/mypanel.vcf.gz` has one haploid column per assembly, which is what kMate
needs, and it is already reduced to top-level bubbles.

Put the jobstore on **node-local disk**. If your reference contains IUPAC codes, use
**Cactus 3.1.0 or newer**; older versions stop with `Non-ACGTN character`.

---

## Step 2: make it biallelic

Graph VCFs are multi-allelic; kMate needs one ALT per record.

> ### Do not use `bcftools norm -m -any`
> On a pangenome graph, splitting records by realignment **scatters carriers across
> shifted positions and silently drops them** where ALT paths converge on the same
> variant. Measured on our panel: up to **~99% carrier loss** at one site. Nothing
> errors; you get wrong frequencies.

Use symbolic-ID propagation instead, which matches variants by identity rather than by
alignment. Give the graph's VCF and GFA:

```bash
kmate decompose \
    --genotyped-vcf out/mypanel.vcf.gz \
    --gfa           out/mypanel.gfa.gz \
    --out panel.vcf.gz
```

The output is sorted and indexed. No `--haploidize`: assembly genotypes are already
haploid.

This is the slow step: it reads the whole GFA (Chr1 of a 135-assembly *Arabidopsis*
graph: ~30 min, ~20 GB). Subset the VCF to one chromosome first (`bcftools view -r Chr1`)
and run chromosomes in parallel; the GFA stays whole.

### Acknowledgement

The decomposition method is not kMate's. It is **HPRC symbolic-ID propagation**, run
via `annotate_vcf.py` and `convert-to-biallelic.py` from
[PanGenie](https://github.com/eblerjana/pangenie) (Jana Ebler, MIT). kMate bundles both
unmodified, so there is nothing extra to install.

If the decomposition matters to your results, cite
Ebler et al. (2022) *Nature Genetics* 54:518–525 and
Liao et al. (2023) *Nature* 617:312–324.

---

## Step 3: build the k-mer index

```bash
kmate build-index --vcf out/mypanel.vcf.gz --ref REF.fa --out index/ours -k 31 --haploid
```

This takes the **graph** VCF from Step 1, not the decomposed one: the index is built per
graph bubble, and a decomposed VCF splits bubbles into overlapping records, which
`build-index` rejects.

Produces `index/ours_<CHR>_kmers.tsv.gz`: the k-mers that identify each bubble.

---

## Step 4: build the two matrices

```bash
# founder x k-mer
kmate build-kmer-pa \
    --kmers index/ours_Chr1_kmers.tsv.gz \
    --vcf panel.vcf.gz --ref REF.fa --chrom Chr1 \
    --out kmer_pa/kmer_pa_Chr1 \
    --treat-missing-as-n --filter-production --min-ac 1 --invariant-margin 1

# founder x variant
kmate build-var-pa --vcf panel.vcf.gz --chrom Chr1 --out var_pa/var_pa_Chr1
```

### `--min-ac` depends on how many founder haplotypes you have

This drops uninformative k-mers. The right value is **not** the same for every panel:

| founder haplotypes | use | why |
|---|---|---|
| many (≈100+) | `--min-ac 2` | a k-mer in only one haplotype is usually noise |
| few (≲20) | `--min-ac 1` | haplotype-private k-mers are your *most* informative ones |

On an 8-founder panel, `--min-ac 2` threw away **~46% of all k-mers**. The default is 2,
so **lower it for a small panel**.

---

## Other starting points

Assemblies of every founder (Step 1) are the case kMate is designed for. The routes
below work, with caveats.

### Some founders genotyped from short reads

If some founders have no assembly, they can be genotyped on the graph (e.g. with
PanGenie) and added to the panel. This is what the GrENE-Net 231-founder panel does:
78 assemblies plus 153 short-read founders.

**Avoid mixing sources if you can.** An assembly resolves far more k-mers than a
short-read genotype of the same founder, so the two kinds of column are not equally
complete. On the GrENE-Net panel this skewed the mixture toward assembled founders
until kMate gained its per-founder normalization; mixed panels need more checking
than all-assembly ones.

Genotyped VCFs lack the graph's variant IDs and are diploid, so decompose them against
an annotation of the graph's own VCF (make it once with `annotate_vcf.py`, bundled in
`kmate/_vendor/pangenie/`), and make the genotypes haploid:

```bash
kmate decompose \
    --annotated-vcf      annotated_multiallelic.vcf.gz \
    --biallelic-catalog  annotated_biallelic.vcf.gz \
    --genotyped-vcf      your_genotyped.vcf.gz \
    --haploidize --het missing \
    --out genotyped_panel.vcf.gz
```

then merge with the assembly side. `--het` depends on what the founders are; see
[Haplotypes and windows](Haplotypes-and-windows):

| your founders | use | what happens |
|---|---|---|
| **inbred lines**: *Arabidopsis* accessions, MAGIC/RIL founders, NAM parents | `--het missing` | `0/0`→`0`, `1/1`→`1`, heterozygous→missing |
| **outbred and phased** | `--het split` | each sample becomes two haplotype columns, `sample.h1` and `sample.h2` |

For inbred founders a heterozygous call is usually an error, so marking it missing lets
kMate skip that founder at that variant instead of guessing. kMate reports the rate and
warns if more than 10% of calls are heterozygous, which usually means your founders are
**not** inbred.

`--het split` needs **phased** genotypes; an unphased heterozygote is an error, because
guessing the phase would invent haplotypes. It is also the route if Cactus wrote two
haplotypes of one individual as a single diploid sample (it does so for names like
`plantA.1`, `plantA.2`).

### An existing pangenome graph

Reduce its VCF to top-level bubbles, then decompose with `--gfa` as in Step 2, using the
graph's GFA:

```bash
vcfbub -l 0 -a 100000 --input graph.vcf.gz | bgzip > panel.vcf.gz
bcftools index -t panel.vcf.gz
```

### A phased multi-sample callset

This works, but you only get what the callset contains: usually SNPs and short
indels, no structural variants.

---

## Check the panel before using it

```python
import numpy as np, scipy.sparse as sp
K = sp.load_npz("kmer_pa/kmer_pa_Chr1.kmer_pa.npz")
V = sp.load_npz("var_pa/var_pa_Chr1.var_pa.npz")
C = sp.load_npz("var_pa/var_pa_Chr1.var_called.npz")

assert K.shape[0] == V.shape[0] == C.shape[0]   # same haplotypes, same order
print("haplotypes:", K.shape[0], "k-mers:", K.shape[1], "variants:", V.shape[1])
print("variants with no called haplotype:", int((np.asarray(C.sum(0)).ravel() == 0).sum()))
```

- The **founder-haplotype count must match** across all three files.
- **Variants with no called haplotype** have undefined frequency; there should be none.
- A haplotype with far fewer k-mers than the rest is usually a broken sample column.
