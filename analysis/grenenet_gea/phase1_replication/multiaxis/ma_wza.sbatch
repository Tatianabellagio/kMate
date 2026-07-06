#!/bin/bash
#SBATCH --job-name=ma_wza
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --array=0-39
#SBATCH --output=/global/scratch/users/tbellg/kmate/logs/ma_wza_%A_%a.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/logs/ma_wza_%A_%a.out
# Multi-axis reblock + WZA (deg-2 canonical, caps snp=1000/nonsnp=350, no floor hack).
# One array task per (axis, class): reblock kendall/lfmm/binomial(quasi) onto clq0.9,
# then WZA each -> multiaxis/wza/wza_{model}_{cls}_gen9_{axis}_deg2.csv
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
source analysis/grenenet_gea/phase1_replication/multiaxis/axes.sh
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PR=analysis/grenenet_gea/phase1_replication
MA=results/grenenet_gea/phase1_replication/multiaxis
i=${SLURM_ARRAY_TASK_ID}; AXIS=$(axis_for $i); C=$(cls_for $i)
CAP=1000; [ "$C" = "nonsnp" ] && CAP=350
echo "== host $(hostname) | ma_wza $C $AXIS (cap=$CAP) =="

echo "-- reblock 3 models onto clq0.9 --"
$PY -u $PR/multiaxis/reblock_multiaxis.py --axis "$AXIS" --cls "$C" --r2 0.9

echo "-- WZA deg2 per model --"
for M in kendall lfmm binomial; do
  IN=$MA/wza_in/${M}_${C}_gen9_${AXIS}.csv
  [ -f "$IN" ] || { echo "  SKIP $M (no $IN)"; continue; }
  $PY -u $PR/run_wza.py --model "$M" --class "$C" --gen 9 --climate "$AXIS" \
      --indir "$MA/wza_in" --out "$MA/wza" --regime deg2 --sample-snps $CAP
done
echo "== ALLDONE_MA_WZA $C $AXIS =="
