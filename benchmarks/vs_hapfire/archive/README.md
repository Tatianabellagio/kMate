# vs_hapfire/archive/

Superseded kMate-vs-hapFIRE material. Nothing here is current.

| dir | what it is |
|---|---|
| `panel_mismatch_prefit_2026-07-15/` | **The important one.** All `hapfire_vs_kmate_*` results from before the panel-mismatch fix: kMate ran on arch3 while hapFIRE was handed the old greneNet SNP-only VCF, and reads came from arch3 — so hapFIRE was solving against a reference its reads did not match (2.46% of reads at NM>10). Any "hapFIRE underperforms kMate on accession recovery" claim from these runs is invalid. Full writeup: [`../PANEL_MISMATCH_BUG.md`](../PANEL_MISMATCH_BUG.md). |
| `single_run_speed_2026-06-03/` | the original single-condition speed run (n231, Chr1, cov50) — superseded by the full N × depth × seed sweep in `../results/`. |
| `diagnostic_tests_2026-07/` | one-off probes: `genomewide_test`, `oldpanel_test`, `phase4_test`, `r1_test`. Kept for provenance of the setup decisions, not for numbers. |
| `plot_consolidation_2026-07-21/` | figure/scorer versions replaced when the plots were consolidated on 2026-07-21. |

The live comparison is the parent dir: speed, founder-`h` accuracy, and
`../ecotype_resolution/` — all on **tool-native panels**, which is the basis the
panel-mismatch bug and the 2026-08-25 shared-panel decision both point to.
