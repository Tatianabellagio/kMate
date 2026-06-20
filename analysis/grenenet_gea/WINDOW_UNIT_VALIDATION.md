# Window-mode block unit — production run, validation, and panel cleanup (2026-06-20)

The selection-test unit is now the **LD-defined block** (dynld K500 unit), and the cohort
has been re-run in **window mode** so every sample carries a per-unit haplotype-frequency
vector `h` (the natural selection-unit state) alongside the projected per-variant AF.

## 1. Production window-mode run
- Runner: `grenenet/run_site_array_perchrom.sh` with `BLOCK_MODE=window`,
  `BLOCKS_DIR=results/grenenet_gea/blocks_mcf90` (units `chr{N}_units_dynld_K500.tsv`,
  22,939 units, 72% k-mer-covered), `MIN_KMERS=50`. Output: `results/grenenet_kmate_window/`.
- 2,168 / 2,168 samples. Per sample: `*_ChrN.h_blocks_per_chrom.npz` (per-unit h) +
  `*_ChrN.tsv` / `*.tsv` (projected per-variant AF). ~1.6 TB.

## 2. Window vs global AF (does the unit-local h change AF?)
`compare_window_vs_global_af.py` + `compare_win_vs_glob.sbatch` + `aggregate_window_vs_global.py`
→ `results/grenenet_gea/window_vs_global/` (SUMMARY.txt, plot).
- 18.4B variant-comparisons, all 2168 samples. **median Pearson r 0.995**, mean |ΔAF| 0.0077.
- ~48% of variants ~identical (desert/global-fallback); **17% shift >0.01, 1.1% >0.10** —
  window redistributes AF within the covered LD units (the haplotype-resolved refinement).

## 3. Panel homogenization (removed the mixed-panel footgun)
Sites 4 & 54 (175 samples) + the 8 SEEDMIX gen0 reps had their **global** TSVs on the OLD
10.33M panel; the rest are on the 8.49M segregating-only panel. Because the filter only
changed V_pa (projection), NOT K_pa (h estimation), correcting = an exact **row-subset**.
- `correct_oldpanel_sample.py` + `correct_oldpanel.sbatch`: line-wise subset via
  `old2new_mask.npy`, originals archived to `*_oldpanel_archive/`. Verified: per-chrom mask
  slices reproduce new-panel positions exactly; kept rows byte-identical; dropped rows
  monomorphic. **All raw TSVs (cohort arch3 + SEEDMIX) are now 8.49M.**
- The old→new mask handling code was then DELETED from `build_af_store.py`, `build_p0.py`,
  `build_two_stage_pooled.py`, `_build_support_nb.py` (they now assert `len == n_full`).

## 4. Block breakage with block-based h (founder vs evolved PC1-VE)
`analysis/grenenet_gea/gen9_window/` (build_gather_index → extract_sample → merge_pools →
founder_vs_evolved_dynld) → notebook `notebooks/block_breakage_window_vs_global.ipynb`.
Recomputes the founder-vs-evolved-PC1-VE breakage diagnostic on the **dynld units**, with the
evolved gen9 PC1-VE from window h vs global AF (same 355 pools + same records).
- 22,879 units: breakage (drop>0.2) **6.3% (window) vs 6.7% (global)** — essentially equal;
  concentrated in desert units (window=global by fallback). Block-based h neither manufactures
  nor hides breakage. NB: absolute evolved VE (~0.50) is lower than the original clq0.9-block
  plot because dynld units are coarser (median 32 var vs 7), not because of the AF mode.

## Next
GEA selection testing on the per-unit `h` (Pipeline B v2, `PIPELINE_B_POOLED_MODEL.md`):
pool plots → site freq → weighted within-site slope → IV climate regression w/ permutation null.
