# arch3/chr1/archive — concluded one-off diagnostics

Superseded exploration and validation one-offs from the Chr1 arch3 build. They
are **not** part of the production build chain (A1–A5, D1) and not the kept
open-loop validation (`../jobA7_compare_vs_hapfire.sh`). Preserved for
archaeology; not maintained, paths are run-specific.

| file | what it was |
|---|---|
| `jobA6_project_h_seedmix.sh` | re-projected a pre-computed SEEDMIX_S1 h-vector through the new var_pa; superseded by `src/per_sample_per_chrom.py` (one pass). |
| `jobB1_seedmix_trim_clumpify.sh` | read-side clumpify dedup for SEEDMIX FASTQs; belongs to read preprocessing (`docs/PIPELINE_FASTQ_PREPROCESSING.md`), not var_pa. |
| `jobC1_outlier_mechanism_decomp.sh` | one-off decomposition of |Δ|>0.10 outliers (encoding vs genuine); conclusion folded into `docs/INVESTIGATION_CN_VAR_DECOMPOSITION.md`. |
| `jobD2_atomized_project_and_compare.sh` | "did atomization help?" experiment (answer: yes). |
| `jobD3_atomized_outlier_deepdive.sh` | atomized-outlier characterization (centromere/F_MISSING strata). |
| `jobE1_af_truth_vs_estimate.sh` + `af_truth_vs_estimate_arch3.py` | 4-panel scatter; its "recipe truth" panel is `uniform_h @ var_pa` (panel-intrinsic, **not** open-loop truth). |
| `jobE2v2_h_vs_recipe_diagnostic.sh`, `jobE2_h_vs_recipe_diagnostic.sh` | closed-loop h-vs-recipe deviation diagnostic (E2 is the slow predecessor of E2v2). Not validation. |
| `jobE3_execute_notebook.sh` | nbconvert wrapper for a notebook; broken path. |
| `jobF1_filt2_project_and_compare.sh` | k-mer-filter sweep (filt2 vs mixedloose); decision closed 2026-05-27 (`docs/METHODS_TRIED_AND_RESULTS.md`). |
