# Archived: Beagle imputation pipeline workdir

Moved here on 2026-05-08 from `imputation/work/` and `imputation/work_merged/`.

**Canonical reasoning lives at the project level — read these instead:**
- `METHODS_TRIED.md` § 5 — "HARD REJECT" entry with the headline numbers
- `RESULTS_LOG.md` 2026-05-03 entry — per-class concordance table (75-LOO test)
- `RESULTS_LOG.md` 2026-05-06 entry — literature corroboration + lock-in decision
- `pangenie_genotyping/data/merged/README_GOLDEN_STANDARD.md` — production VCF wired to use the pre-imputation catalog

This directory is a candidate for deletion once disk is needed. Nothing here
is referenced by production code paths; the rejection rationale survives in
the project-level docs above.

## What's in this directory

```
work/        Beagle pipeline intermediates: subset extracts (grene_{80,151,all}),
             per-step VCFs (cactus_svs_renamed, target_151, target_svs_missing,
             ref_80, imputed_151), cn matrices (cn_231/), Chr4 test_fix workdir,
             per-ecotype LOO Beagle runs (loo/<ecotype>/ × 11)

work_merged/ Per-chrom Beagle outputs (imputed_chr{1..5}.vcf.gz), the merged
             biallelic + multiallelic Beagle output (founders_231_imputed*.vcf.gz),
             and the pre-Beagle prep VCF (founders_231_dipl_split.vcf.gz)
```

The scripts that produced these files stayed in `imputation/` (small, document
the approach). One reference to a moved path remains in
`imputation/test_fix_chr4.sh` — historical, not production.

## How to restore

```bash
mv archive/imputation_beagle/work        imputation/work
mv archive/imputation_beagle/work_merged imputation/work_merged
```
