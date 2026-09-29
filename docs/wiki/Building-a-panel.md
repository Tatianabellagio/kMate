# Building a panel

Do this **once** per founder set. Afterwards you reuse the matrices for every pool.

```
your founders ──▶ panel VCF ──▶ k-mer index ──▶ kmer_pa
                       └──────────────────────▶ var_pa + var_called + meta
```

All four steps are `kmate` commands. Run them **per chromosome**.

---

## Step 1 — get a panel VCF

The VCF must meet the requirements in [Input files](Input-files#input-variants-the-panel-vcf).
Pick whichever starting point matches what you have.

### A. You have founder assemblies

Build a pangenome graph. This gives the full variant spectrum — SNPs, indels **and**
structural variants — with real genotypes.

```bash
cactus-pangenome <jobstore> panel.seqfile \
    --outDir out --outName mypanel \
    --reference REFNAME --haplo --vcf --gfa --gbz --giraffe
```

`panel.seqfile` is `<name><TAB><path>`, with the **reference first**. The output
`out/mypanel.vcf.gz` is already reduced to top-level bubbles.

Two practical points: put the jobstore on **node-local disk**, not shared storage; and
if your reference contains IUPAC codes, use **Cactus 3.1.0 or newer** — older versions
stop with `Non-ACGTN character`.

### B. You have a pangenome graph VCF already

Reduce it to top-level bubbles:

```bash
vcfbub -l 0 -a 100000 --input graph.vcf.gz | bgzip > panel.vcf.gz
bcftools index -t panel.vcf.gz
```

### C. You have a phased multi-sample callset

This works, but you only get what the callset contains — usually SNPs and short
indels, no structural variants.

---

## Step 2 — make it biallelic and haploid

Graph VCFs are usually multi-allelic. kMate needs one ALT per record and one allele
per genotype.

> ### Do not use `bcftools norm -m -any`
> On a pangenome graph, splitting records by realignment **scatters carriers across
> shifted positions and silently drops them** where ALT paths converge on the same
> variant. Measured on our panel: up to **~99% carrier loss** at one site. Nothing
> errors — you just get wrong frequencies.

Use symbolic-ID propagation instead, which matches variants by identity rather than by
alignment:

```bash
kmate decompose \
    --annotated-vcf      annotated_multiallelic.vcf.gz \
    --biallelic-catalog  annotated_biallelic.vcf.gz \
    --genotyped-vcf      your_genotyped.vcf.gz \
    --convert-to-biallelic /path/to/convert-to-biallelic.py \
    --haploidize --het missing \
    --out panel.vcf.gz
```

### Choosing `--het`

kMate founders are haplotypes, so diploid genotypes must be reduced to one allele.
How depends on **what your founders are**:

| your founders | use | what happens |
|---|---|---|
| **inbred lines** — *Arabidopsis* accessions, MAGIC/RIL founders, NAM parents | `--het missing` | `0/0`→`0`, `1/1`→`1`, heterozygous→missing |
| **outbred and phased** — e.g. HPRC assemblies | `--het split` | each sample becomes two haplotype columns, `sample.h1` and `sample.h2` |

For inbred founders a heterozygous call is usually an error, so marking it missing lets
kMate skip that founder at that variant instead of guessing. kMate reports the rate and
warns if more than 10% of calls are heterozygous — that usually means your founders are
**not** inbred and you want `--het split` instead.

`--het split` needs **phased** genotypes. An unphased heterozygote is an error, not a
coin flip, because guessing the phase would invent haplotypes.

### Attribution

`kmate decompose` mostly **runs tools written by other people**, which are not bundled —
you obtain them separately. Run `kmate decompose --citation` to print this.
**If you publish results that used it, cite them:**

| tool | from | cite |
|---|---|---|
| `annotate_vcf.py` | Human Pangenome Reference Consortium | Liao et al. (2023) Nature 617:312–324 |
| `convert-to-biallelic.py` | [eblerjana/pangenie-tools](https://github.com/eblerjana/pangenie-tools) | Ebler et al. (2022) Nature Genetics 54:518–525 |
| `bcftools` | samtools project | Danecek et al. (2021) GigaScience 10:giab008 |

Only the `INFO/ID` transfer step is kMate's own.

---

## Step 3 — build the k-mer index

```bash
kmate build-index --vcf panel.vcf.gz --ref REF.fa --out index/ours -k 31 --haploid
```

Produces `index/ours_<CHR>_kmers.tsv.gz`: the k-mers that identify each bubble.

---

## Step 4 — build the two matrices

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

### `--min-ac` depends on how many founders you have

This drops uninformative k-mers. The right value is **not** the same for every panel:

| founders | use | why |
|---|---|---|
| many (≈100+) | `--min-ac 2` | a k-mer in only one founder is usually noise |
| few (≲20) | `--min-ac 1` | founder-private k-mers are your *most* informative ones |

On an 8-founder panel, `--min-ac 2` threw away **~46% of all k-mers**. The default is 2,
so **lower it for a small panel**.

---

## Check the panel before using it

```python
import numpy as np, scipy.sparse as sp
K = sp.load_npz("kmer_pa/kmer_pa_Chr1.kmer_pa.npz")
V = sp.load_npz("var_pa/var_pa_Chr1.var_pa.npz")
C = sp.load_npz("var_pa/var_pa_Chr1.var_called.npz")

assert K.shape[0] == V.shape[0] == C.shape[0]      # same founders, same order
print("founders:", K.shape[0], "k-mers:", K.shape[1], "variants:", V.shape[1])
print("variants with no called founder:", int((np.asarray(C.sum(0)).ravel() == 0).sum()))
```

- The **founder count must match** across all three files.
- **Variants with no called founder** have undefined frequency; there should be none.
- A founder with far fewer k-mers than the rest is usually a broken sample column.
