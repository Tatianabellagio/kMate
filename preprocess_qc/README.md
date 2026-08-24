# preprocess_qc/ — FASTQ preprocessing QC and genotyped-panel concordance

QC of the sequencing cohort as it moves through preprocessing (trim → dedup → align →
genotype), plus concordance of the genotyped panel against GrENE-Net SNP calls.

Written 2026-08-24 during the root-level cleanup.

## ⚠ This directory is a *remainder* — read this before reusing anything

Roughly half of what used to be here was archived on 2026-07-10 as **pre-arch3 and stale**.
`archive/preprocess_qc/` holds the LD-decay / SV-SNP tagging analysis and the panel
composition descriptive stats; both ran on the old merged panel
(`founders_231_chr.vcf.gz`), not on production arch3. Their current replacements are:

| Archived analysis | Current replacement |
|---|---|
| SNP/SV LD-decay + SV-SNP tagging (`notebooks/ld_*.ipynb`, `scripts/{compute_ld,aggregate_ld,compute_per_sv_max_r2}.py`, `output/ld/`) | `analysis/grenenet_gea/` SV-SNP tagging work → `analysis/grenenet_gea/sv_snp_ld/` |
| Panel composition stats (`notebooks/production_vcf_stats.ipynb`, `scripts/merged_vcf_stats*`) | `analysis/panel_qc/panel_stats/PANEL_STATS.md` + `scripts/panel_stats_for_paper.py` |

**What remains here is the non-LD, non-composition part** — duplicate rates, trim survival,
stage-by-stage retention, genotyped-sample QC, and GrENE-Net concordance. That part was
never declared stale.

## Scripts

**Preprocessing-stage QC**

| Script | What it does |
|---|---|
| `qc_one.sh` | Per-sample QC for one library (SLURM array task); writes to `output/qc_results/` (git-ignored, regeneratable). |
| `build_qc_index.py` | Assembles per-task results into `qc_index_{main,loo,redo}.tsv`. |
| `parse_trim_stats.py` | Trimmomatic survival → `trim_stats_{main,loo}.tsv`. |
| `parse_dup_rates.py` | Clumpify duplicate rates → `dup_rates_{main,loo}.tsv`. |
| `build_stage_aggregate.py` | Stage-by-stage retention per ecotype → `stage_aggregate.tsv`. |
| `audit_presence.py` | Which samples are actually present at each stage → `presence_{main,loo}.tsv`. |

**Genotyped-panel QC and concordance**

| Script | What it does |
|---|---|
| `qc_genotyped_one.py` / `qc_genotyped_array.sh` | Per-sample QC of the genotyped panel → `genotyped_qc/`, aggregated to `genotyped_qc_main_aggregate.tsv`. |
| `grenenet_snp_concordance.py` / `grenenet_concord_array.sh` | Concordance vs GrENE-Net SNP calls → `grenenet_concordance/`, aggregated to `grenenet_concordance_aggregate.tsv` and `concordance_combined.tsv`. |

## Outputs

Small aggregate TSVs at `output/*.tsv` and the four retention/coverage figures in
`output/plots_stages/` are **tracked**. The bulky per-task and per-sample trees
(`output/qc_results/`, `output/genotyped_qc/`, `output/grenenet_concordance/`) and `logs/`
are git-ignored and regeneratable — see this directory's own `.gitignore` plus the
`preprocess_qc/output/…` rules in the repo-root `.gitignore`.

## Notebooks

`notebooks/preprocess_qc.ipynb` (cohort QC overview) and
`notebooks/preprocess_qc_stages.ipynb` (stage-by-stage retention). Executed copies
(`*.executed.ipynb`) are git-ignored.

Narrative write-up: `docs/PIPELINE_FASTQ_PREPROCESSING.md`.
