# blocks — what a test unit *is*

Defines the analysis unit that r1/r2/r3 aggregate over. This is a small
**package**, not a pile of scripts: `recompute_blocks` and
`eval_block_coherence` each have 8 importers, and `block_cluster_pc1ve`
(`cluster_founders`) and `block_haplotype_counts` (`hap_counts`) are imported by
three scripts apiece. Do not move or retire a file here without checking
`from <name> import` first — several read as one-off Chr1 prototypes but are
libraries the production path depends on.

## Map construction

| script | what it builds |
|---|---|
| `recompute_blocks.py` | the faithful HapFM/BigLD recompute on our all-class founder `var_pa` |
| `blocks_tiling.py` | the **gap-free "tiling"** assignment (hapFIRE's own rule) — the production partition |
| `dynamic_ld_blocks.py` | LD-guided growth: relax r² only where blocks are too small |
| `make_coarse_blocks.py` | greedy adjacent merge to coarsen a map |
| `build_hap_membership.py` | founder → haplotype membership (the projection matrix M_b) |

sbatch arms: `blocks_recompute{,_mcf90,_p80}`, `blocks_clq05_mcf90`,
`clq_sweep_mcf90`, `cluster_mcf90`, `dynld`, `hap_membership{,_clq50,_clq90,_clq90nosv}`,
`diag_genomewide`, `unit_dist`, `install_bigld`. These are **parameter arms, not
generations** — `blocks_mcf90` is the locked production map.

> ⚠ **The tiling `merge_small_blocks` step has a known failure mode on sparse
> classes.** For SVs (~27k genome-wide) most tiling blocks hold 0–1 SVs, so
> merging chains empty blocks backward into the last SV-dense one. An SV-class
> block's genomic footprint can reach tens of kb — far wider than its ~1–2 kb LD
> island — which silently attributes a lead SV to distant, unlinked genes. This
> is what produced the CARK8/9 mis-attribution. **Never treat an SV block
> footprint as its LD island**; check the lead variant's own position
> (`genes/dissection/`).

## Coherence diagnostics — does a block move as one unit?

`blockcoherence_data.py` (PC1-VE per block across CLQcut levels) ·
`delta_p_coherence.py` (direction-concordance metrics on Δp, not just static LD) ·
`eval_block_coherence.py` (shared evaluator, imported by 8) ·
`coherence_vs_floor.py` (coherence cost of coarsening) ·
`block_founder_vs_evolved.py` + `verify_founder_ve.py` (founder vs evolved PC1-VE,
the recombination diagnostic) · `block_unit_frontier.py` (the two-decision frontier)

## Characterization

`block_haplotype_counts.py` (how multi-haplotype are the blocks, on founders) ·
`block_cluster_pc1ve.py` (HapFM x-means clustering of low-VE blocks; exports
`cluster_founders`) · `block_kmer_coverage.py` (panel k-mer count per block) ·
`block_panel_support_tag.py` (per-variant tagging + panel support) ·
`block_missing_sensitivity.py` (did imputation distort the clusters?) ·
`unit_distributions.py` + `plot_unit_distributions.py` ·
`build_sv_landscape.py` + `plot_sv_landscape.py` (per-block SV landscape and the
two confounds any enrichment must survive) ·
`_check_clq90_blocks.py` (disjointness + **inter-block gap fraction** — the loss
that motivated the tiling partition)

Notebooks: `block_coherence_clqcut.ipynb`, `blocks_units_decision.ipynb`,
`haploblock_r20{10,20}_eps0_validation.ipynb`,
`block_breakage_window_vs_global.ipynb` (retired window-mode subsystem, kept as
record). Standing decisions: `WINDOW_UNIT_VALIDATION.md`.
