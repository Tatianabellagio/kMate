# Input files

To estimate frequencies for one pool, kMate needs **two things**: the pooled reads,
and a **panel** describing your founders.

```
pooled reads  ─┐
               ├──▶  kmate run  ──▶  allele frequencies
panel matrices ┘
```

---

## Input reads

- **Paired FASTQ** from one pool (`R1.fq.gz`, `R2.fq.gz`). A single-end file works too.
- **Adapter- and quality-trimmed.**
- **Deduplicated, if your library is PCR-based.** PCR duplicates inflate k-mer counts.
  PCR-free libraries need no deduplication.
- One pool per run. Run kMate once per sample.

Coverage: a few × is enough for usable estimates. Below ~1× the result is dominated
by sampling noise — check coverage and drop failed libraries **before** estimating,
not after.

---

## The panel

The panel is built **once** per founder-haplotype set and reused for every pool. It is
four files per chromosome:

| file | what it holds |
|---|---|
| `kmer_pa_<CHR>.kmer_pa.npz` | which founder haplotypes carry each k-mer — the evidence for the mixture |
| `var_pa_<CHR>.var_pa.npz` | which founder haplotypes carry the ALT at each variant — the projection target |
| `var_pa_<CHR>.var_called.npz` | which founder haplotypes have a genotype call at each variant |
| `var_pa_<CHR>.meta.npz` | the variants themselves: chromosome, position, REF, ALT |

`kmer_pa` also has a small `.meta.npz` beside it.

**You build these yourself from your founders** — see
[Building a panel](Building-a-panel). It is a one-time job per founder set; every pool
you sequence afterwards reuses the same matrices.

> The matrices are large (tens of GB for a few hundred founders) and are not stored in
> the git repository. The 231-founder *Arabidopsis thaliana* panel used by GrENE-Net is
> available on request.

---

## Input variants (the panel VCF)

Only relevant if you are **building** a panel. Both matrices are built from one VCF,
and it must satisfy all of the following.

- **Multi-sample** — one column per **founder haplotype**. The sample order becomes the
  haplotype order and is shared by both matrices.
- **Haploid** — exactly one allele per genotype (`0`, `1`, `.`), not `0/1`.
  A panel column *is* a haplotype — see [Haplotypes and windows](Haplotypes-and-windows).
  The builders stop with an error on diploid input.
- **Biallelic** — one ALT per record. The builders stop with an error otherwise.
- **Sequence-resolved** — real REF/ALT sequences. Symbolic alleles like `<DEL>` cannot
  be turned into k-mers.
- **Non-overlapping** — one variant per locus per haplotype. For graph VCFs this means
  top-level bubbles only (what `vcfbub` produces).
- **bgzip-compressed and indexed** — `.vcf.gz` plus `.tbi` or `.csi`.
- **Segregating only**, recommended — drop variants carried by nobody or by everybody;
  they cost memory and add nothing.

[Building a panel](Building-a-panel) shows how to get from a raw VCF to one that
satisfies these.

---

## Input reference

- A **FASTA** matching your VCF's coordinates and contig names, indexed (`.fai`).
- **One reference for every step.** IUPAC ambiguity codes are fine — kMate converts
  any non-ACGTN base to `N` internally, the same thing a pangenome graph build does.
