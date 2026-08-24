#!/bin/bash
#SBATCH --job-name=man_kendall
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --array=0-19
#SBATCH --output=/global/scratch/users/tbellg/kmate/logs/man_kendall_%A_%a.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/logs/man_kendall_%A_%a.out
# Pooled non-SNP (sv U smallindel) Kendall-tau, gen9, one array task per axis.
# For the pooled snp-vs-nonsnp before/after-WZA table (faithful nonsnp on the
# current-AF nonsnp class matrix; the 3-class run stays the primary product).
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
source analysis/grenenet_gea/phase1_replication/multiaxis/axes.sh
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
AXIS=${AXES[$SLURM_ARRAY_TASK_ID]}
OUT=analysis/grenenet_gea/phase1_replication/results/multiaxis/kendall
echo "== host $(hostname) | kendall nonsnp gen9 $AXIS =="
$PY -u analysis/grenenet_gea/phase1_replication/run_kendall.py \
    --class nonsnp --gen 9 --climate "$AXIS" --threads "${SLURM_CPUS_PER_TASK:-8}" --out "$OUT"
echo "== ALLDONE_MAN_KENDALL nonsnp $AXIS =="
