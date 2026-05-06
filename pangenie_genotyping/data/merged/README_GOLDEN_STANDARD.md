# GOLDEN-STANDARD production VCF

```
founders_231_chr.vcf.gz   ← THIS IS THE PRODUCTION DELIVERABLE
founders_231_chr.vcf.gz.tbi
```

- **231 founder samples**: 80 cactus-assembly genotypes (long-read truth) + 151 PanGenie short-read calls
- **5,214,959 records**: SVs + small indels + SNPs in one catalog
- **Mixed ploidy**: cactus 80 are haploid (`0`, `1`, `.`); PanGenie 151 are diploid (`0/0`, `0/1`, `1/1`)
- **Multi-allelic preserved** (cactus pangenome bubbles, not biallelic-decomposed)
- **Cactus-side haploid `.` preserved** — biologically meaningful ("this assembly's path doesn't traverse this bubble"), not "missing by quality"

## Why no imputation

Tested Beagle 5.5 imputation (xwu's pattern) on this VCF in 2026-05-03. Result:
- SNP concordance: −0.5pp (noise)
- **small_sv concordance: −29.9pp**
- **medium_sv concordance: −24.7pp**
- 0 of 75 LOO ecotypes improved post-imputation

Confirmed by literature (cattle pangenome paper, French dairy SV imputation
study, PanGenie original): SV imputation under-performs direct PanGenie
genotyping at coverage ≥10×. Production workflows in similar projects keep
PanGenie SV calls unimputed.

See `RESULTS_LOG.md` 2026-05-06 entry for full literature backing.

## How downstream tools should consume this

`build_cn_var.py` already handles the mixed-ploidy + cactus-`.` correctly:
- Per-record: `has_alt = any(a > 0 for a in gt)` → `cn=1` for any carrier (haploid or diploid), `cn=0` for hom_ref or `.`/`./.`
- For inbred *A. thaliana*, this carrier-status interpretation is biologically appropriate (per 2026-05-03 het-rate analysis showing PanGenie het is mostly artifact).

If a downstream tool needs imputed SNPs (e.g., hapFIRE-style LD analyses),
produce a SEPARATE SNP-only Beagle-imputed VCF — do NOT co-impute SVs.

## Related files (NOT used in production)

- `imputation/work_merged/founders_231_imputed.vcf.gz` — Beagle-imputed (biallelic-decomposed). Reference only.
- `imputation/work_merged/founders_231_imputed_multiallelic.vcf.gz` — re-merged to multi-allelic. Reference only.
- `imputation/work_merged/founders_231_dipl_split.vcf.gz` — pre-Beagle prep step. Intermediate.

## Stats locations

- Per-record / per-sample / per-size-class TSVs: `preprocess_qc/output/merged_stats/founders_231_*.tsv*`
- Visual notebook: `preprocess_qc/notebooks/production_vcf_stats.ipynb`
