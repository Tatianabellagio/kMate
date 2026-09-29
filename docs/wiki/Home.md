# kMate

kMate estimates **allele frequencies for SNPs, short indels and large structural
variants** from pooled sequencing of a population descended from a known set of
**founder haplotypes**.

It is **alignment-free**: reads are never mapped. Working in k-mer space is why
structural variants are handled exactly like SNPs.

---

## How it works, in two steps

### Step 1: put your founder haplotypes into k-mer space

Turn your founder haplotypes' pangenome into **two matrices**, once per founder set and
reused for every pool.

A *founder haplotype* is one continuous sequence, one allele per variant. An inbred line
is one haplotype; a phased outbred individual is two. See
[Haplotypes and windows](Haplotypes-and-windows).

```
                                    ┌─ K_pa : founder × k-mer
   founders' pangenome (a VCF) ─────┤
                                    └─ V_pa : founder × variant
```

| matrix | what it says | what it is for |
|---|---|---|
| **K_pa** | for each k-mer, which founder haplotypes contain it | the **evidence**: what the reads are compared against |
| **V_pa** | for each variant, which founder haplotypes carry the ALT allele | the **translation**: turns haplotype proportions into allele frequencies |

A third file, `var_called`, records where each founder haplotype actually had a call, so
that missing data is excluded rather than silently counted as reference.

→ **[Building a panel](Building-a-panel)**

### Step 2: run a pool against them

Count the k-mers in the pooled reads and solve for the mixture that produced them: each
count is modelled as Poisson with mean `λ · Σ h_f · K_pa[f,k]`, and the founder-haplotype
mixture **`h`** is fitted by expectation-maximisation.

Projecting `h` through **V_pa** gives each variant's frequency: the summed proportion of
haplotypes carrying it. One pass covers every SNP, indel and SV.

```
   pooled reads ──▶ k-mer counts ──▶ EM ──▶ h (haplotype mixture) ──▶ × V_pa ──▶ allele frequencies
```

→ **[Running kMate](Running-kMate)**

---

## Start here

| | |
|---|---|
| **[Installation](Installation)** | install kMate and check it works |
| **[Input files](Input-files)** | what kMate needs, and what each file must satisfy |
| **[Haplotypes and windows](Haplotypes-and-windows)** | what a founder haplotype is, and whether to run per chromosome or per window |
| **[Building a panel](Building-a-panel)** | **step 1**: build the matrices from your founder haplotypes |
| **[Running kMate](Running-kMate)** | **step 2**: estimate frequencies for a pool |
| **[Output](Output)** | the result table, and how to read it |
| **[Command reference](Command-reference)** | every option, with its default |
| **[Troubleshooting](Troubleshooting)** | common errors and what they mean |

---

## The short version

```bash
# install
mamba create -n kmate -c conda-forge -c bioconda python numpy scipy pysam kmer-jellyfish samtools
mamba activate kmate
pip install kmate
kmate selftest                      # must print PASS
```

**Step 1: build the panel from your founders' VCF** (once):

```bash
kmate build-index    --vcf panel.vcf.gz --ref REF.fa --out index/ours -k 31 --haploid

kmate build-kmer-pa  --kmers index/ours_Chr1_kmers.tsv.gz \
                     --vcf panel.vcf.gz --ref REF.fa --chrom Chr1 \
                     --out kmer_pa/kmer_pa_Chr1 \
                     --treat-missing-as-n --filter-production

kmate build-var-pa   --vcf panel.vcf.gz --chrom Chr1 --out var_pa/var_pa_Chr1
```

**Step 2: run each pool against them** (per sample):

```bash
kmate run \
    --kmer-pa-prefix kmer_pa/kmer_pa \
    --var-pa-prefix  var_pa/var_pa \
    --reads R1.fq.gz R2.fq.gz \
    --sample MYPOOL --out MYPOOL_Chr1.tsv \
    --chroms Chr1 --unit chrom
```

Step 1 requires a haploid, biallelic, sequence-resolved VCF.
[Building a panel](Building-a-panel) covers how to get one, including from raw assemblies.

---

## Source, issues, citation

Code and the detailed technical docs: <https://github.com/Tatianabellagio/kMate>

kMate was developed for the [GrENE-Net](https://www.science.org/doi/10.1126/science.adz0777)
outdoor evolution experiment in *Arabidopsis thaliana*. It is not yet published as a
standalone method. If you want to use it, or collaborate, get in touch:
**Tatiana Bellagio** (tatianabellagio@gmail.com).

The panel-building path builds on tools from HPRC and eblerjana; see
[Building a panel](Building-a-panel#acknowledgement).
