# ⛔ Benchmark AF-scoring rules — READ BEFORE comparing tools

Two mistakes were made on 2026-06-20 that produced a fake "kMate is worse at allele
frequency than hapFIRE" result. Both are easy to repeat. Don't.

## RULE 1 — Never join AF estimates by `(chrom, pos)` alone

`recomb_truth.tsv.gz` and the kMate `_global.tsv` carry **no REF/ALT bases** — only
`ref_len`/`alt_len`. ~45k Chr1 positions are **multiallelic** (>1 biallelic SNP record at
the same POS). A `(chrom,pos)`-only join (or `groupby(pos).first()`) pairs the **wrong
allele** between a tool's estimate and the truth → garbage R².

- Proof: single-founder pool `cov10_n1_g0` → buggy join gave kMate AF R²=**0.937**; correct
  join gives **1.0000** (physics-required — one founder, AF = its genotype).

**Correct way:** restrict ALL tools to the **shared SNP panel**
(`work/shared_snps_<panel>_Chr1.vcf.gz`). Multiallelic positions were already dropped when
that panel was built, so within it `(chrom,pos)` IS unique and unambiguous — which is the
only reason hapFIRE/vg (which output `chrom,pos,freq` with no allele) can be joined safely.
Use **`accuracy_vs_competitors/scripts/score_snp_fair.py`** — it does this. The old
`score_competitor_snp.py` and `score_competitors.py` are **WRONG** (position-only join).

## RULE 2 — `recomb_truth.truth_af` is MAR; it is a closed loop that flatters kMate

`compute_recomb_truth.py` builds `truth_af` as `Σw·varpa / Σw·called` — the **MAR**
convention (divide by *called* founders), explicitly "to match how cactus_em projects AF."
So kMate's estimand == the truth's convention. Scoring **any all-records / low-info** subset
against `truth_af` structurally **favors kMate** and **sinks read-based hapFIRE/vg** (which
estimate *physical* AF — the reads were built missing→REF). hapFIRE R²=−0.39 at info<0.5 is
*partly estimand mismatch, not error*.

**Headline only the convention-free basis: fully-called SNPs (`info ≥ 0.99`),** where
MAR ≡ physical. There kMate ≈ hapFIRE (e.g. n80_g0: 0.9946 vs 0.9978). Same caveat for SVs:
the honest number is fully-called R²=0.992, **not** the all-SV 0.756 (that tail is the MAR
confound). Never headline MAR / all-records.

## Bottom line for the paper
SNP AF is a **tie** with hapFIRE (not a win, not a loss). kMate's real, defensible wins are
**speed, alignment-free, SVs, and ecotype/`h` resolution** (the last is projection-free, so
it has neither of the above confounds). See memory `snp-af-scoring-bug-and-mar-confound`.
