# Archived: window-unit / haploblock-frequency (hapfreq) diagnostics — RETIRED 2026-07-08

## What this is
Scripts from the retired **window-unit** and **haploblock-frequency (hapfreq)**
branch of the GrENE-net GEA. They estimated per-*window* founder frequencies and
projected them into per-haplotype-cluster frequencies, plus a family of
window-vs-global / cross-chromosome agreement diagnostics.

## Why it's retired
1. Production kMate now runs `--unit chrom` (GLOBAL mode; see
   `../../GLOBAL_MODE_DECISION.md`). The `--unit chrom` cohort produces **no
   per-block `h_blocks`** files.
2. These scripts read the **window stores**
   `results/grenenet_kmate_window` and `results/grenenet_kmate_window_seedmix`
   (per-sample `*_Chr{N}.h_blocks_per_chrom.npz`). **Those directories have been
   deleted from disk**, so every script here is broken/unrunnable as-is.
3. The window-vs-chrom question is settled: at r²=0.1 with eps=0 the haploblocks
   collapse to ~231 genome-wide (≈ chromosome-wise), so the window/haploblock
   reframe is "framing only" — it does not change the AF/`h` estimates. Earlier
   agreement reads (r≈0.99–0.995) were consistent with this but predate the
   full-panel `Kf_w` normalization fix and were partly circular (see
   `../../WINDOW_UNIT_VALIDATION.md`, `../../GLOBAL_MODE_DECISION.md`).

## Contents
- `cross_chrom_agreement.py`, `cross_chrom_by_generation.py` — cohort cross-chromosome
  `h`-agreement diagnostics (read the deleted window store).
- `block_residual_consistency.py`, `block_vs_global_h.py` — per-site block-vs-global
  `h` consistency diagnostics.
- `compare_analyses_site.py` — per-site window-vs-global comparison.
- `compare_win_vs_glob.sbatch` — driver (also referenced the old
  `results/grenenet_kmate_arch3` store, which is NOT retired and stays on disk).
- `gen9_window/` — the full window/gen9 "breakage" diagnostic script subtree
  (clustering, extraction, seedmix, haploblock/dynld founder-vs-evolved, notebook
  builders). Self-contained: no code outside this subtree imported it.
  **NOTE:** the associated DATA directory `results/grenenet_gea/gen9_window/`
  (kept for GEA) was **NOT** moved — only the scripts are archived here.

## NOT archived (deliberately left in place — see the handoff report)
The core hapfreq/Pipeline-B chain — `build_hapfreq_matrix.py`,
`build_hapfreq_p0_seedmix.py`, `build_hap_trajectories.py`, `build_hap_gea.py`,
`build_hap_wza.py`, and its drivers `hapfreq.sbatch` / `gea_clq90_pipelineB.sbatch`,
plus `phase1_replication/build_hap_lastgen_matrix.py` — was left in place. Although
these also read the deleted window store (`build_hapfreq_matrix.py` reads
`results/grenenet_kmate_window` even in `--h-source global` mode), they generate
`results/grenenet_gea/hapfreq_clq90/pipelineB_varlen/*` which is still consumed by
the **live** SV-enrichment analysis (`sv_adaptive/sv_enrichment_gea.py`). Retiring
this chain requires a maintainer decision and should be done as a unit.

Date archived: 2026-07-08.
