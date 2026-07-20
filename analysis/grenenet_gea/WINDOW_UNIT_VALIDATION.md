# Block/window units — definitions, production run, validation (2026-06-19/20)

> **Status note (added 2026-07-06).** Window-mode kMate as an **AF estimator** is **superseded** by
> `GLOBAL_MODE_DECISION.md` (2026-07-01): evolved allele frequencies are estimated in GLOBAL mode.
> The **LD-defined blocks / dynld units below are retained** — not for AF estimation, but as the
> **GEA test units** (the unit selection acts on; median n_eff ≈ 2.46 haplotypes/block). The
> window-mode cohort run and its window-vs-global AF comparison (§2–§4) stand as the validation that
> window-local `h` does not materially change AF (median r 0.995) — i.e. GLOBAL loses little.
> This file merges `BLOCKS_HANDOFF.md` (block definitions + benchmark handoff); original in git `5cefcfa`.
>
> **Archival note (added 2026-07-08).** The window-unit / haploblock-frequency (hapfreq)
> *diagnostic* machinery this file describes has been **archived** to
> `archive/window_hapfreq_retired/` (cross-chrom agreement, block-vs-global-h, the
> `gen9_window/` script subtree). The window stores `results/grenenet_kmate_window[_seedmix]`
> those scripts read have been **deleted from disk**; the `--unit chrom` production cohort
> emits no per-block `h_blocks`. The r≈0.995 window-vs-global agreement numbers below
> **predate the full-panel `Kf_w` normalization fix** and should be read as historical.
> Window-vs-chrom is settled: haploblocks (r²=0.1, eps=0) collapse to ~231 ≈ chromosome-wise,
> so the reframe is "framing only." See `GLOBAL_MODE_DECISION.md`.

The selection-test unit is the **LD-defined block** (dynld K500 unit); the cohort was also
re-run in **window mode** so every sample carries a per-unit haplotype-frequency vector `h`
alongside the projected per-variant AF.

## 0. Block / unit definitions (the unit map)

**The units** (k-mer-covered, LD-grown; = kMate h-window candidate AND selection unit):
- genome-wide: `analysis/grenenet_gea/blocks_mcf90/final_units_dynld_K500.tsv`; per-chrom
  `chr{N}_units_dynld_K500.tsv`. Cols: `chrom start_pos end_pos n_variants panel_kmers covered`.
  **22,939 units; 72% COVERED** (panel_kmers≥500 → local-fit; 16,403 units), 28% desert (global
  fallback). Median 32 var / 2.1 kb / 926 panel k-mers. Built by `dynamic_ld_blocks.py` (grow CLQ0.9
  blocks along the LD gradient until ≥500 panel k-mers). Code-audited + validated (2026-06-19,
  3-agent adversarial review; k-mer tagging byte-identical to `block_em.assign_kmers_to_blocks`).
- Base CLQ0.9 blocks (before dynld growth): `chr{N}_clq0.9_blocks_clq0.9.tsv`, **58,376 blocks**
  genome-wide, median 223 bp / 7 variants (many too thin for h-estimation → hence the dynld K500 grow).
  Built by `recompute_blocks.py` (HapFM `CompleteLDPartition` corr=0.2 + `BigLD` gpart) via
  `blocks_recompute_mcf90.sbatch` (`--corr 0.2 --clqcut 0.9 --min-called-frac 0.9`), on the all-class
  panel `panel/arch3/chr{N}/var_pa_231_arch3_chr{N}.*`.
- Haplotype units within blocks: each n_eff≤2 block = 1 unit (block AF); each n_eff>2 block split
  into HapFM-xmeans clusters keeping PC1-VE≥0.7 (`block_cluster_pc1ve.py`, needs the `np.warnings`
  shim; `block_unit_frontier.py`). Registry `final_units_ve07.csv` (66,032 units).
- Coarsened floor maps (greedy-merge to a min-variant floor; `make_coarse_blocks.py`):
  `blocks_mcf90/coarse/floor{8,15,25,40}.tsv`. Env gotchas: BigLD R needs
  `export LD_LIBRARY_PATH=/usr/lib64`; xmeans needs the numpy-warnings shim; run heavy jobs via sbatch.
  Characterization: `notebooks/blocks_units_decision.ipynb`.

## 1. Production window-mode run
- Runner: `grenenet/run_site_array_perchrom.sh` with `BLOCK_MODE=window`,
  `BLOCKS_DIR=analysis/grenenet_gea/blocks_mcf90` (units `chr{N}_units_dynld_K500.tsv`,
  22,939 units, 72% k-mer-covered), `MIN_KMERS=50`. Output: `results/grenenet_kmate_window/`.
- 2,168 / 2,168 samples. Per sample: `*_ChrN.h_blocks_per_chrom.npz` (per-unit h) +
  `*_ChrN.tsv` / `*.tsv` (projected per-variant AF). ~1.6 TB.

## 2. Window vs global AF (does the unit-local h change AF?)
`compare_window_vs_global_af.py` + `compare_win_vs_glob.sbatch` + `aggregate_window_vs_global.py`
→ `analysis/grenenet_gea/archive/window_hapfreq_retired/window_vs_global/` (SUMMARY.txt, plot).
- 18.4B variant-comparisons, all 2168 samples. **median Pearson r 0.995**, mean |ΔAF| 0.0077.
- ~48% of variants ~identical (desert/global-fallback); **17% shift >0.01, 1.1% >0.10** —
  window redistributes AF within the covered LD units (the haplotype-resolved refinement).

## 3. Panel homogenization (removed the mixed-panel footgun)
Sites 4 & 54 (175 samples) + the 8 SEEDMIX gen0 reps had their **global** TSVs on the OLD
10.33M panel; the rest are on the 8.49M segregating-only panel. Because the filter only
changed V_pa (projection), NOT K_pa (h estimation), correcting = an exact **row-subset**.
- `correct_oldpanel_sample.py` + `correct_oldpanel.sbatch`: line-wise subset via
  `old2new_mask.npy`, originals archived to `archive/*_oldpanel_archive/` (moved from
  `results/*_oldpanel_archive/` 2026-07-10, see `archive/README.md`). Verified: per-chrom mask
  slices reproduce new-panel positions exactly; kept rows byte-identical; dropped rows
  monomorphic. **All raw TSVs (cohort arch3 + SEEDMIX) are now 8.49M.**
- The old→new mask handling code was then DELETED from `build_af_store.py`, `build_p0.py`,
  `build_two_stage_pooled.py`, `_build_support_nb.py` (they now assert `len == n_full`).

## 4. Block breakage with block-based h (founder vs evolved PC1-VE)
`analysis/grenenet_gea/archive/window_hapfreq_retired/gen9_window/` (build_gather_index → extract_sample → merge_pools →
founder_vs_evolved_dynld) → notebook `notebooks/block_breakage_window_vs_global.ipynb`.
Recomputes the founder-vs-evolved-PC1-VE breakage diagnostic on the **dynld units**, with the
evolved gen9 PC1-VE from window h vs global AF (same 355 pools + same records).
- 22,879 units: breakage (drop>0.2) **6.3% (window) vs 6.7% (global)** — essentially equal;
  concentrated in desert units (window=global by fallback). Block-based h neither manufactures
  nor hides breakage. NB: absolute evolved VE (~0.50) is lower than the original clq0.9-block
  plot because dynld units are coarser (median 32 var vs 7), not because of the AF mode.

## Next

Not window mode. Pipeline B v2 (`PIPELINE_B_POOLED_MODEL.md`) depended on a per-unit
window-mode `h`, which is dead: given ~97% selfing + ~3 generations, a window-mode rerun
is not worth doing (`GLOBAL_MODE_DECISION.md`). GEA selection testing proceeds on
`--unit chrom` (global-mode) AF only.
