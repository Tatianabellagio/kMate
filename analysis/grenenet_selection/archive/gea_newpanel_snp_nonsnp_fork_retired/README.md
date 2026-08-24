# New-panel GEA — SNP vs non-SNP (raw LFMM + WZA), clq0.9 haploblocks

The kMate arch3 panel adds **indels and SVs** on top of the SNP-only phase-1 paper.
This folder asks, on that new panel: **run a raw climate-GEA (LFMM) separately on
SNPs vs pooled non-SNPs (indel+SV) and see how similar the results are — where do
the hits land, and do non-SNPs reveal signal SNPs miss?**

Split off from `../phase1_replication/` on 2026-07-02 (that folder stays as the
phase-1 CAM5 replication). This folder holds only the SNP-vs-nonSNP comparison.

## Key choices
- **Model:** raw **LFMM ridge, K=16**, `lfmm_test(calibrate='gif')`, env = bio1.
  K=16 confirmed by GIF sweeps on both classes (GIF plateaus ~2.1–2.5, never →1;
  the GIF-calibrated null is flat at every K → K=16 = phase-1 convention is fine).
- **Classes:** `snp` (1.99M records) vs `nonsnp` = indel+SV pooled (693k records),
  MAF ≥ 0.05, gen9 last-gen pools (355 `site_gen_plot`, flower-weighted Δp).
- **Blocks:** **clq0.9 haploblocks** (`blocks_clq09.py`, from
  `results/.../blocks_recompute/chr*_clq0.9_blocks_clq0.9.tsv`) — BigLD CLQ-cut 0.9,
  ~82k blocks genome-wide (median 188 bp / 7 variants), far finer than the phase-1
  hapFIRE LD blocks (16.7k). We **do not** use the phase-1 blocks here.
- **Aggregation:** canonical WZA deg7-cap2000 (`../wza_script.py`) on the clq0.9
  blocks, per class — shown alongside the raw per-record scan.

## Files
- `blocks_clq09.py` — clq0.9 block loader + `assign_clq09_blocks(chrom, pos)`.
- `run_wza_clq09.py` (+`.sbatch`) — reassign clq0.9 blocks to the raw LFMM p-values,
  run WZA per class → `results/.../gea_newpanel/{lfmm_*_clq09.csv, wza_*_clq09.csv}`.
- `_build_snp_vs_nonsnp_nb.py` → `notebooks/snp_vs_nonsnp_lfmm.ipynb` — the deliverable:
  raw-LFMM Manhattan **and** WZA Manhattan (SNP vs non-SNP), raw top hits, and
  gene-annotated class-specific blocks. Built with kMate env, executed in `basic`.
- `run_nonsnp_ksweep.sbatch`, `run_snp_vs_nonsnp_k16.sbatch` — how the non-SNP class
  matrix / LFMM input / K-sweep / K=16 run were produced (provenance).
- **Other two phase-1 tests (raw, no WZA):** `run_{kendall,binomial}_snp_vs_nonsnp.sbatch`
  drive `../phase1_replication/run_{kendall,binomial}.py --class {snp,nonsnp} --gen 9
  --climate bio1` → `results/.../gea_newpanel/{kendall,binomial}/{test}_{cls}_gen9_bio1.csv`.
  `_build_binom_kendall_nb.py` → `notebooks/binom_kendall_snp_vs_nonsnp.ipynb` — raw
  per-record Manhattans (binomial + Kendall, SNP vs non-SNP), gene-annotated top hits,
  and concordance (binomial↔kendall within class; SNP↔non-SNP LD-block peak per test).
  (`nonsnp` was added to the `--class` choices of both phase-1 runners.)

## Upstream inputs (shared, live in `../phase1_replication/`)
- `build_class_matrices.py` (with the added `nonsnp` class) → per-class pool matrices.
- `build_lfmm_input.py` (with `nonsnp`) + `run_lfmm_lastgen.R` → raw LFMM p per record
  (`results/.../phase1_replication/lfmm/lfmm_{snp,nonsnp}_gen9_bio1.csv`).

## Outputs → `analysis/grenenet_gea/gea_newpanel/results/`
- `lfmm_{cls}_gen9_bio1_clq09.csv` — per-record LFMM p with clq0.9 block.
- `wza_{cls}_clq09.csv` — per-block WZA (`gene`=block, `Z_pVal`, chrom, pos).
- `snp_vs_nonsnp/` — Manhattan PNGs, top-record + class-specific block tables (annotated).
