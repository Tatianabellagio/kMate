Event-level GEA scan (sv + nonsnp), 2026-09-30: MAF >= 0.05 + call rate >= 0.50, SVs as Truvari events.

For analysis read the consolidated production tables instead:
  r2_gea_nonsnp/phase1_replication/results/multiaxis/tile_cr50/   (all 4 classes, one schema, shared tiling blocks)
The lfmm_event_*.csv.gz here hold the same p-values (block ids differ: these are hapFIRE blocks, tile_cr50 uses the clq0.9 tiling).

Unique to this folder (read these here):
  sv_gen9_units.csv, nonsnp_gen9_units.csv   unit definitions after the call-rate floor
  sv_gen9_dp.npy, nonsnp_gen9_dp.npy         event-level frequency-change matrices the scan used
  lambda_by_axis.csv, event_sv_hits.csv, event_direction*.csv   derived summaries
The *_units.csv one level up (genes/convergence/event_gea/) are the PRE-floor unit sets, used on purpose by
annotate_event_units.py and _sv_context_decomposition.py as the "before the floor" set.
