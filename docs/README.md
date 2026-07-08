# docs/

Analyses, investigations, and methods writeups for the kMate project. Entry-point
docs (`README.md`, `HANDOFF.md`, `ALGORITHM.md`, `BACKGROUND.md`, `SAVIO_HPC.md`)
live at the project root.

| File | Purpose |
|---|---|
| `PIPELINE_STATE.md` | **Single source of truth** for what is production vs in-evaluation vs deprecated (panel, projection, filters). Read before making pipeline decisions. |
| `METHODS_TRIED_AND_RESULTS.md` | What was tried and how it scored — k-mer-filter sweep, EM weighting (ω_k=1/m_b vs uniform), projection variants. Complements `RESULTS_LOG.md`. |
| `FOUNDER_NORMALIZATION_FIX.md` | The per-founder M-step normalization fix (`--normalize per_founder`, full-panel Kf_w) that stopped the global-mode founder-h collapse; supersedes ω_k=1/m_b for the selfing production path. |
| `EM_UNIT_CHOICE_AND_NONIDENTIFIABILITY.md` | Why `--unit chrom` is the production estimator (and `--unit ld` collapses in low-diversity blocks); EM non-identifiability and the no-prior decision. |
| `RERUN_AFTER_FIX.md` | Re-run checklist for regenerating all downstream outputs after the per_founder + Kf_w fix (now landed in `results/grenenet_gea/rerun_kfw_hb`). |
| `INVESTIGATION_CN_VAR_DECOMPOSITION.md` | Why production switched to the arch decomposition (annotate_vcf + convert-to-biallelic) instead of `bcftools norm -m -any` for var_pa. |
| `PIPELINE_FASTQ_PREPROCESSING.md` | Read-side FASTQ preprocessing pipeline (trim/dedup/QC) used upstream of genotyping. |
| `MISSINGNESS_231PANEL.md` | F_MISSING characterization on the production 231-founder panel and its effect on AF estimates. |
| `SIMULATIONS_METHODS.md` | Methods-ready description of the pool-seq VISOR simulation framework (regime matrix, parameters, citations). |
| `EXTERNAL_DEPENDENCIES.md` | What the gitignored `external/` dir must contain (upstream URLs + pinned commits) and which production step uses each — so a fresh clone can reconstruct it. |
| `RESULTS_LOG.md` | Chronological numerical record (large file). Useful historical reference; not authoritative for current state — see `PIPELINE_STATE.md`. |
