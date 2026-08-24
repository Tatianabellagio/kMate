#!/bin/bash
#SBATCH --job-name=ma_quasibinom
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=16
#SBATCH --mem=110G
#SBATCH --time=8:00:00
#SBATCH --array=0-59
#SBATCH --output=/global/scratch/users/tbellg/kmate/logs/ma_quasibinom_%A_%a.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/logs/ma_quasibinom_%A_%a.out
# Multi-axis per-variant quasi-binomial + K=16 latent factors (the production
# binomial model), gen9. One array task per (axis, class).
# Output -> multiaxis/quasibinom/quasibinom_lf16_{cls}_gen9_{axis}.csv (pval_quasi col)
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
source analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/multiaxis/axes.sh
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python   # statsmodels in basic
PR=analysis/grenenet_selection/phase1_replication   # production quasibinom now lives here
OUT=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/multiaxis/quasibinom
i=${SLURM_ARRAY_TASK_ID}; AXIS=$(axis_for $i); C=$(cls_for $i)
echo "== host $(hostname) | quasibinom $C K16 gen9 $AXIS =="
$PY -u $PR/run_quasibinom.py --class "$C" --gen 9 --climate "$AXIS" \
    --k 16 --threads "${SLURM_CPUS_PER_TASK:-16}" --out "$OUT"
echo "== ALLDONE_MA_QUASIBINOM $C $AXIS =="
