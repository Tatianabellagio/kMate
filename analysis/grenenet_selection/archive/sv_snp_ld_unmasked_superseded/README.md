# Retired: SV-only, unmasked SV-SNP tagging (r2) analysis

Archived 2026-07-10 (`build_sv_snp_ld.py` only — its notebook, `sv_snp_ld_tagging.ipynb`,
was gone from disk at the time). The notebook was later regenerated from the same old
buggy `.npz` outputs via `_build_sv_snp_ld_nb.py` and sat live in `notebooks/` presenting
the retracted numbers below as current findings; both are now archived alongside it
(2026-07-20). These computed max founder-genotype r2
between each SV and its best nearby SNP (arch3 panel vs GrENE-Net short-read),
but had two bugs found while reviewing the plots:

1. **Missingness silently coded as REF.** Only `var_pa` (ALT-carrier) was loaded,
   never `var_called` (the called mask) — a missing genotype and a true REF call
   were indistinguishable, for both the SV and the candidate SNP.
2. **No MAC floor.** 53.5% of SVs are singletons (MAC=1) — a singleton SV could
   register a spurious r2=1 "tag" against any other marker private to the same
   one founder, a small-N coincidence, not real LD.

`SUMMARY.md` (in `output/`) documents the old headline numbers (85.8% / 47.4% at
r2>0.8) — **do not cite these**, they predate both fixes.

**Superseded by:**
- `analysis/grenenet_gea/build_tagging_masked.py` — corrected computation
  (proper pairwise-complete masked r2 + configurable MAC floor), extended to
  indels, output in `analysis/grenenet_gea/sv_snp_ld_v2/`.
- `analysis/grenenet_gea/notebooks/sv_indel_tagging_masked.ipynb` — corrected,
  expanded write-up (SV + indel, panel + short-read, plus missingness-filter
  and MAF-filter sensitivity plots).

The shared inputs this old script also used (`greneNet_final_v1.1.vcf.gz`,
`shortread_geno_Chr{1-5}.npz`) were **not** archived — they still live in
`analysis/grenenet_gea/sv_snp_ld/` and are reused by the corrected script.
