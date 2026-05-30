# GOLDEN-STANDARD production VCF + haploid sibling (SUPERSEDED — v1, historical)

> **Superseded.** This describes the **v1** mixed-ploidy catalog
> (`founders_231_chr.*`, ~5.2M records, built via `bcftools norm -m -any` +
> `haploidize_merged_vcf.sh` — now in `scripts/archive/`). The current production
> panel is **`data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz`** (het-masked,
> AC=0-cleaned; see `../../README.md`), and per-record matrices are built by the
> arch decomposition (`panel/arch3/`), **not** `norm -m -any`. Kept for the durable
> rationale below — the **no-imputation** decision and the **carrier-status (het→.)
> collapse** still hold for v3qc_v3 — but treat the file names, record counts, and
> build steps here as historical.

This directory holds the production 231-founder catalog in two on-disk forms:

```
founders_231_chr.vcf.gz           ← PRODUCTION DELIVERABLE (mixed-ploidy)
founders_231_chr.haploid.vcf.gz   ← haploid sibling (uniform {0, 1, .} per cell)
```

Both files are kept on disk. The mixed-ploidy version remains the on-disk source of truth (cactus haploid + PanGenie diploid as written by the cactus pangenome and the PanGenie genotyper respectively). The haploid sibling is for downstream consumers that prefer a uniform `{0, 1, .}` per-cell encoding.

## Contents (both files)

- **231 founder samples**: 80 cactus-assembly genotypes (long-read truth) + 151 PanGenie short-read calls
- **5,214,959 records**: SVs + small indels + SNPs in one catalog
- **Multi-allelic preserved** (cactus pangenome bubbles, not biallelic-decomposed)
- **Cactus-side haploid `.` preserved** — biologically meaningful ("this assembly's path doesn't traverse this bubble"), not "missing by quality"

The mixed-ploidy file has cactus 80 as haploid (`0`, `1`, `.`) and PanGenie 151 as diploid (`0/0`, `0/1`, `1/1`); the haploid file collapses both sides to `{0, 1, .}` (see "Carrier-status collapse" below).

## Why no imputation

Tested Beagle 5.5 imputation (xwu's pattern) on this VCF in 2026-05-03. Result:
- SNP concordance: −0.5pp (noise)
- **small_sv concordance: −29.9pp**
- **medium_sv concordance: −24.7pp**
- 0 of 75 LOO ecotypes improved post-imputation

Confirmed by literature (cattle pangenome paper, French dairy SV imputation study, PanGenie original): SV imputation under-performs direct PanGenie genotyping at coverage ≥10×. Production workflows in similar projects keep PanGenie SV calls unimputed.

See `docs/RESULTS_LOG.md` 2026-05-06 entry for full literature backing.

## How downstream tools should consume this

`build_cn_var.py` handles the mixed-ploidy + cactus-`.` correctly:

- Per-record: `has_alt = any(a > 0 for a in gt)` → `cn=1` for any carrier (haploid or diploid), `cn=0` for hom_ref or `.`/`./.`
- For inbred *A. thaliana*, this carrier-status interpretation is biologically appropriate (per 2026-05-03 het-rate analysis showing PanGenie het is mostly artifact — see "Carrier-status collapse rationale" below).

If a downstream tool needs imputed SNPs (e.g., hapFIRE-style LD analyses), produce a **separate SNP-only Beagle-imputed VCF** — do NOT co-impute SVs.

---

## Haploid sibling — how built

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

Script: `panel/pangenie_genotyping/scripts/haploidize_merged_vcf.sh`.

## Carrier-status collapse rationale

**Het rate is small and dominated by PanGenie call noise on a known subset of samples — not real biology.** Direct count on the production VCF (Chr1:1-2,000,000, all 151 PanGenie samples × 81,084 records = 12.2 M calls, 2026-05-08):

| GT category | count | % of all calls |
|---|---|---|
| `0/0` hom_ref | 11,093,556 | **90.6 %** |
| `1/1` hom_alt | 929,591 | **7.6 %** |
| `0/1` het | 96,214 | **0.79 %** |
| `1/2`, `0/2`, … multi | 123,344 | 1.0 % |
| `./.` missing | 979 | 0.008 % |

Of alt-carrying calls (n = 1,149,149), **92 % are already 1/1** and only 8.4 % are het. Per-sample distribution (n = 151):

| percentile | het % |
|---|---|
| P25 | 0.45 % |
| P50 (median) | 0.62 % |
| P75 | 0.87 % |
| P95 | 1.95 % |
| max | 5.06 % (ecotype 9977) |

144 of 151 samples are below 2 % het. The seven samples ≥ 2 % het are exactly the Cao 2011 GAII short-read libraries (9977, 9985, 10013, 9634, 9748, 9743, 9941) flagged for noisy PanGenie calls in the 2026-05-03 RESULTS_LOG entry — their "het" cells correlate r ≈ 0.96 with disagreement against independent GrENE-Net SNP truth, confirming they are PanGenie ambiguity, not biological heterozygosity.

For inbred *A. thaliana*, the true biological het rate is ≤ 0.5 %. Collapsing PanGenie diploid genotypes to haploid carrier-status reclassifies < 1 % of total calls and ≈ 10 % of alt-carrying calls as 1 instead of 0/1 — and that reclassification is biologically defensible (inbred carriers ARE essentially homozygous; the `0/1` cells are mostly call noise, not real heterozygotes).

## What downstream consumers gain from the haploid file

- LD / r² computations no longer need a mixed-ploidy → carrier-dose mapping trick — every cell is `{0, 1, .}`.
- Multi-allelic decomposition is precomputed (one pass, persisted on disk) instead of repeated by every consumer.
- Per-record encoding is uniform between the cactus 80 and the PanGenie 151 sides — no mixed-ploidy footnote in any pipeline.
- Tools that don't natively handle mixed ploidy can consume the haploid VCF directly.

What is preserved in the haploid file:
- Cactus-side haploid `.` is preserved as `.` (not promoted to `./.` or imputed) — same biological semantics.
- PanGenie-side `./.` is preserved as `.`.
- Every record from the source VCF that produced ≥1 biallelic record is represented.
- No imputation. The 2026-05-06 RESULTS_LOG decision still holds.

## Related files (NOT used in production)

- `panel/imputation/work_merged/founders_231_imputed.vcf.gz` — Beagle-imputed (biallelic-decomposed). Reference only.
- `panel/imputation/work_merged/founders_231_imputed_multiallelic.vcf.gz` — re-merged to multi-allelic. Reference only.
- `panel/imputation/work_merged/founders_231_dipl_split.vcf.gz` — pre-Beagle prep step. Intermediate.

## Note for hapFIRE consumers

Neither file is directly consumable by hapFIRE. hapFIRE requires phased diploid GTs (regex `[0-9]\|[0-9]`) and one biallelic record per `(chrom, pos)` (HARP per-base likelihood is biallelic). Both files here are unphased / multi-allelic-atomized.

For the hapFIRE methods-comparison column on the v3 panel, the conversion pipeline is:

1. `sims/visor_freqk/scripts/build_chr1_panel_v3.sh` — Chr1 biallelic SNPs; chrom rename `Chr1`→`1`; collapse to phased homozygous diploid carrier-status (`0|0`/`1|1`); rebuild fixed-window block-index.
2. `sims/visor_freqk/scripts/dedup_panel_v3.sh` — keep one biallelic record per `(chrom, pos)`; rebuild block-index against the dedup'd panel.

See `sims/visor_freqk/chr1_only_panel_v3/README.md` for the full processing write-up.

## Stats locations

- Per-record / per-sample / per-size-class TSVs: `preprocess_qc/output/merged_stats/founders_231_*.tsv*`
- Visual notebook: `preprocess_qc/notebooks/production_vcf_stats.ipynb`

## Audit trail

- 2026-05-08 — Haploid sibling built from `founders_231_chr.vcf.gz` (5,214,959 records, 231 samples). After `bcftools norm -m -any`, ALT count differs because multi-allelic records split.
- Het analysis underlying the carrier-status decision: `docs/RESULTS_LOG.md` 2026-05-03 entry ("PanGenie het rate is mostly artifact, not biology") + 2026-05-08 confirmation count above.
- Imputation rejection: `docs/RESULTS_LOG.md` 2026-05-06 entry (literature lock-in).
