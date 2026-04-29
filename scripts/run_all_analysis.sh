#!/bin/bash
# Run the full simulation analysis once hapFIRE jobs are done.
set -euo pipefail
cd /home/tbellagio/scratch/hapfire_sv

PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python

echo "=== 1. Simulation comparison (rep9, cov50, 1kb) ==="
$PYTHON scripts/plot_simulation.py 2>&1 | tee logs/plot_simulation.log

echo
echo "=== 2. Seedmix comparison (already done; rerun for completeness) ==="
$PYTHON scripts/seedmix_truth_analysis.py 2>&1 | tail -30

echo
echo "=== 3. Quick stats summary ==="
echo "  Outputs in: $(pwd)/results/"
ls -la results/*.png results/*.tsv results/*.csv.gz 2>&1 | head -20
