# Retired: old-panel (arch3 multinomial) scripts

Archived 2026-07-08. These scripts read the pre-fix old-panel kMate output dirs
`results/grenenet_kmate_arch3` (evolved) and `results/seedmix_kmate_arch3`
(seedmix), which were produced under the **global multinomial normalization**
(before the full-panel `Kf_w` / `--normalize per_founder` fix) and have been
**deleted** (~1.6 TB reclaimed). Production AF now lives in
`results/grenenet_gea/rerun_kfw_hb/{evolved,seedmix}` (`--unit chrom` + Kf_w).

- `correct_oldpanel.sbatch` — one-off remap of old-panel per-sample outputs
  (required `af_store/old2new_mask.npy`, which no longer exists). Superseded by
  the full cohort re-run.
- `block_breakage_pilot.py` — Chr1 mosaic/breakage pilot that read old-panel
  `.tsv` AF.
- `crk_block_af_support.py` — phase-1 followup reading old-panel AF.

Do not resurrect against the deleted dirs; repoint to `rerun_kfw_hb` if ever needed.
