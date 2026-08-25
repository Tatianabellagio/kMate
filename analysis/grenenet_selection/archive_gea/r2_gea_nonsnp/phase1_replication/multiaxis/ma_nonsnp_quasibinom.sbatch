#!/bin/bash
#SBATCH --job-name=man_qbinom
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=16
#SBATCH --mem=110G
#SBATCH --time=8:00:00
#SBATCH --array=0-19
#SBATCH --output=/global/scratch/users/tbellg/kmate/logs/man_qbinom_%A_%a.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/logs/man_qbinom_%A_%a.out
# Pooled non-SNP per-variant quasi-binomial + K=16 latent factors, gen9, per axis.
# Output -> multiaxis/quasibinom/quasibinom_lf16_nonsnp_gen9_{axis}.csv (pval_quasi)
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
source analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/multiaxis/axes.sh
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python   # statsmodels in basic
PR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication
OUT=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/multiaxis/quasibinom
AXIS=${AXES[$SLURM_ARRAY_TASK_ID]}
echo "== host $(hostname) | quasibinom nonsnp K16 gen9 $AXIS =="
$PY -u $PR/run_quasibinom.py --class nonsnp --gen 9 --climate "$AXIS" \
    --k 16 --threads "${SLURM_CPUS_PER_TASK:-16}" --out "$OUT"
echo "== ALLDONE_MAN_QUASIBINOM nonsnp $AXIS =="
