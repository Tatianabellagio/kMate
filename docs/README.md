# docs/

Analyses, investigations, and methods writeups for the kMate project. Entry-point
docs (`README.md`, `HANDOFF.md`, `ALGORITHM.md`, `BACKGROUND.md`, `SAVIO_HPC.md`)
live at the project root.

| File | Purpose |
|---|---|
| `PIPELINE_STATE.md` | **Single source of truth** for what is production vs in-evaluation vs deprecated (panel, projection, filters). Read before making pipeline decisions. |
| `METHODS_TRIED_AND_RESULTS.md` | What was tried and how it scored — k-mer-filter sweep, EM weighting (ω_k=1/m_b vs uniform), projection variants. Complements `RESULTS_LOG.md`. |
| `INVESTIGATION_CN_VAR_DECOMPOSITION.md` | Why production switched to the arch decomposition (annotate_vcf + convert-to-biallelic) instead of `bcftools norm -m -any` for var_pa. |
| `PIPELINE_FASTQ_PREPROCESSING.md` | Read-side FASTQ preprocessing pipeline (trim/dedup/QC) used upstream of genotyping. |
| `MISSINGNESS_231PANEL.md` | F_MISSING characterization on the production 231-founder panel and its effect on AF estimates. |
| `SIMULATIONS_METHODS.md` | Methods-ready description of the pool-seq VISOR simulation framework (regime matrix, parameters, citations). |
| `RESULTS_LOG.md` | Chronological numerical record (large file). Useful historical reference; not authoritative for current state — see `PIPELINE_STATE.md`. |
