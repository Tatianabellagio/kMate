# Haploid sibling of the production VCF

```
founders_231_chr.haploid.vcf.gz   ← carrier-status haploid (this file)
founders_231_chr.vcf.gz           ← original mixed-ploidy (preserved, see README_GOLDEN_STANDARD.md)
```

Both files are kept on disk. The haploid version is for downstream consumers
that prefer a uniform `{0, 1, .}` per-cell encoding; the mixed-ploidy version
remains the on-disk source of truth (cactus haploid + PanGenie diploid as
written by the cactus pangenome and PanGenie genotyper respectively).

## How the haploid version is built

```
founders_231_chr.vcf.gz
  │
  ▼  bcftools norm -m -any
multi-allelic records split into biallelic
(cactus integer alleles 33, 51, … and PanGenie 1/2 calls collapse cleanly)
  │
  ▼  awk on FORMAT=GT
   0/0, 0|0, 0    → 0     (hom_ref, cactus haploid 0)
   0/1, 1/0, 1|0,
   0|1, 1/1, 1|1, 1 → 1   (any-alt carrier; cactus haploid 1)
   ./., .|., .    → .     (missing)
  │
  ▼  bgzip + tabix
founders_231_chr.haploid.vcf.gz
```

Script: `pangenie_genotyping/scripts/haploidize_merged_vcf.sh`.

## Why carrier-status haploid is the right collapse for this panel

**Het rate is small and dominated by PanGenie call noise on a known subset of
samples — not real biology.** Direct count on the production VCF
(Chr1:1-2,000,000, all 151 PanGenie samples × 81,084 records = 12.2 M calls,
2026-05-08):

| GT category | count | % of all calls |
|---|---|---|
| `0/0` hom_ref | 11,093,556 | **90.6 %** |
| `1/1` hom_alt | 929,591 | **7.6 %** |
| `0/1` het | 96,214 | **0.79 %** |
| `1/2`, `0/2`, … multi | 123,344 | 1.0 % |
| `./.` missing | 979 | 0.008 % |

Of alt-carrying calls (n = 1,149,149), **92 % are already 1/1** and only 8.4 %
are het. Per-sample distribution (n = 151):

| percentile | het % |
|---|---|
| P25 | 0.45 % |
| P50 (median) | 0.62 % |
| P75 | 0.87 % |
| P95 | 1.95 % |
| max | 5.06 % (ecotype 9977) |

144 of 151 samples are below 2 % het. The seven samples ≥ 2 % het are exactly
the Cao 2011 GAII short-read libraries (9977, 9985, 10013, 9634, 9748, 9743,
9941) flagged for noisy PanGenie calls in the 2026-05-03 RESULTS_LOG entry —
their "het" cells correlate r ≈ 0.96 with disagreement against independent
GrENE-Net SNP truth, confirming they are PanGenie ambiguity, not biological
heterozygosity.

For inbred *A. thaliana*, the true biological het rate is ≤ 0.5 %. Collapsing
PanGenie diploid genotypes to haploid carrier-status reclassifies < 1 % of
total calls and ≈ 10 % of alt-carrying calls as 1 instead of 0/1 — and that
reclassification is biologically defensible (inbred carriers ARE essentially
homozygous; the `0/1` cells are mostly call noise, not real heterozygotes).

## Practical impact

**No information loss for downstream consumers.** `build_cn_var.py` already
applies the carrier-status rule in memory (`has_alt = any(a > 0 for a in
gt)`), so the haploid on-disk file matches what the cn-builder would have
produced anyway. Same goes for the LD computation in
`preprocess_qc/scripts/compute_ld.py` (already maps mixed ploidy → carrier
dose) and for the per-SV r² computation in `compute_per_sv_max_r2.py`.

**What becomes simpler:**

- LD / r² computations no longer need the haploid → diploid `{0, 2}` mapping
  trick — every cell is `{0, 1, .}`.
- Multi-allelic decomposition is precomputed (one pass, persisted on disk)
  instead of repeated by every consumer.
- Per-record encoding is uniform between the cactus 80 and the PanGenie 151
  sides — no mixed-ploidy footnote in any pipeline.
- Tools that don't natively handle mixed ploidy (some VCF readers, some
  population-genetics R packages) can consume the haploid VCF directly.

**What is preserved:**

- Cactus-side haploid `.` is preserved as `.` (not promoted to `./.` or
  imputed) — same biological semantics as in the mixed VCF: "this assembly's
  path didn't traverse this bubble".
- PanGenie-side `./.` is preserved as `.`.
- Every record from the source VCF that produced ≥ 1 biallelic record is
  represented in the haploid VCF.
- No imputation is applied. The Beagle decision in RESULTS_LOG 2026-05-06
  ("don't impute SVs") still holds.

## Audit trail

- 2026-05-08 — Built. Source: `founders_231_chr.vcf.gz` (5,214,959 records,
  231 samples). After `bcftools norm -m -any`, ALT count differs because
  multi-allelic records split.
- Het analysis underlying the decision: `RESULTS_LOG.md` 2026-05-03 entry
  ("PanGenie het rate is mostly artifact, not biology — sticking with
  carrier-status cn") + a fresh count on 2026-05-08 (numbers above) confirming
  the conclusion still holds for the current production VCF.
- Imputation rejection: `RESULTS_LOG.md` 2026-05-06 entry, summarized in
  `README_GOLDEN_STANDARD.md`.
