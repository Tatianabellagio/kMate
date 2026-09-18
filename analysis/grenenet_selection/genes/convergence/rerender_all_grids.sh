#!/bin/bash
# Re-render every candidate figure in the new plot_theme style.
#
# The figures are produced at render time, so a theme change does not touch the ~1,500
# PNGs already on disk -- they have to be redrawn. plot_shortlist.py draws both the locus
# figure and the own-axis garden grid for one candidate table; the reviewed genes are
# spread over three tables, so all three are run. 3 workers each (the node has 4 cores),
# tables sequentially so the workers do not oversubscribe.
#
# ~60-80 s per gene, ~340 genes -> roughly 2-3 hours.
set -u
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
cd "$(dirname "$0")"
mkdir -p results/render_logs
for spec in "screen_top_loci_ownaxis.csv" "screen_gwas_rescreen_v2.csv" "_round1_for_render.csv"; do
  tbl="results/$spec"
  [ -f "$tbl" ] || { echo "skip missing $tbl"; continue; }
  echo "=== $(date +%H:%M) $spec ==="
  for w in 0 1 2; do
    $PY plot_shortlist.py --table "$tbl" --worker $w --nworkers 3 \
        > "results/render_logs/rerender_$(basename $spec .csv)_w$w.log" 2>&1 &
  done
  wait
  echo "=== $(date +%H:%M) $spec done ==="
done
$PY organize_figures.py
echo "=== $(date +%H:%M) ALL DONE, tree re-synced ==="
