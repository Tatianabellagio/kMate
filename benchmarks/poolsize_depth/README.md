# benchmarks/poolsize_depth/

**Question:** how do kMate's accuracy and compute scale with pool size (N founders)
and sequencing depth?

The p231 N × depth sweep — the kMate-side counterpart of the hapFIRE-paper Fig-S13
design. Grid: N ∈ {2,5,20,50,150,231} × depth ∈ {1,10} (+ N=50 × {30,50}) × seeds
42–46, on the arch3 panel.

This dir is **kMate-only**. It is the upstream source for two comparisons that live
elsewhere, which is why it looks small — most of its value is consumed downstream:

| consumer | what it takes from here |
|---|---|
| `../vs_hapfire/` | `results/kmate_speed_table.tsv` (kMate's speed/compute column) and the same N × depth grid definition |
| `../vs_hapfire_ecotype_resolution/` | the matching kMate `h_per_chrom.npz` runs under `../p231/results/kmate_chrom_poolsize_depth/` |

## Scripts

| script | what it does |
|---|---|
| `score_and_plot.py` | the N × depth accuracy sweep → `poolsize_depth_{af,h}_accuracy.png`, `poolsize_depth_table.tsv` |
| `score_snp_vs_nonsnp.py` | kMate AF accuracy split SNP vs non-SNP (indel/SV) on the same grid. Reports the ≥90%-called (`n_called ≥ 208`) filtered values by default — that filtered basis is what `vs_hapfire` consumes as kMate's side |
| `recover_kmate_speed.py` | reconstructs kMate speed/compute for the sweep from existing run logs, **without** re-running the estimator |

Note the estimator outputs themselves live under `../p231/results/kmate_chrom_poolsize_depth/`,
not here; this dir holds only the scored tables and figures.
