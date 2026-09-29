# kMate

kMate estimates **allele frequencies for SNPs, short indels and large structural
variants** from pooled sequencing of a population descended from a known set of
founders.

It works in two steps:

1. From the k-mers in your pooled reads, estimate the **founder mixture** `h` —
   how much of the pool each founder contributes.
2. Project `h` through a founder × variant matrix to get an **allele frequency for
   every variant**, in one pass.

Because the evidence is k-mers rather than aligned bases, SNPs and structural
variants are handled the same way.

---

## Start here

| | |
|---|---|
| **[Installation](Installation)** | install kMate and check it works |
| **[Input files](Input-files)** | what kMate needs, and what each file must satisfy |
| **[Building a panel](Building-a-panel)** | make the matrices from your own founders |
| **[Running kMate](Running-kMate)** | estimate frequencies for a pool |
| **[Output](Output)** | the result table, and how to read it |
| **[Troubleshooting](Troubleshooting)** | common errors and what they mean |

If someone has already given you a panel, you only need
**Installation → Running kMate → Output**.

---

## The short version

```bash
# install
mamba create -n kmate -c conda-forge -c bioconda python numpy scipy pysam kmer-jellyfish samtools
mamba activate kmate
pip install kmate
kmate selftest          # must print PASS

# run one pool against an existing panel
kmate run \
    --kmer-pa-prefix panel/kmer_pa/kmer_pa \
    --var-pa     panel/var_pa/var_pa_Chr1.var_pa.npz \
    --var-called panel/var_pa/var_pa_Chr1.var_called.npz \
    --var-meta   panel/var_pa/var_pa_Chr1.meta.npz \
    --reads R1.fq.gz R2.fq.gz \
    --sample MYPOOL --out MYPOOL_Chr1.tsv \
    --chroms Chr1 --unit chrom
```

---

## Source, issues, citation

Code and the detailed technical docs: <https://github.com/Tatianabellagio/kMate>

kMate was developed for the [GrENE-Net](https://www.science.org/doi/10.1126/science.adz0777)
outdoor evolution experiment in *Arabidopsis thaliana*. It is not yet published as a
standalone method — if you want to use it, or collaborate, get in touch:
**Tatiana Bellagio** (tatianabellagio@gmail.com).

Parts of the panel-building path use tools written by others; see
[Building a panel](Building-a-panel#attribution) for who to cite.
