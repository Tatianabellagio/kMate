#!/bin/bash
#SBATCH --job-name=binom_gea
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/binom_gea_%x_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/binom_gea_%x_%j.out
# Phase-1 binomial-regression GEA, last-gen (gen9), one class per job.
# Usage: sbatch --job-name=binom_<class> run_binomial.sh <class> [gen] [climate]
#   e.g. sbatch --job-name=binom_snp run_binomial.sh snp 9 bio1
set -euo pipefail
CLS=${1:?need class: snp|sv|smallindel}
GEN=${2:-9}
CLIM=${3:-bio1}
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python   # statsmodels lives in basic
PR=analysis/grenenet_gea/phase1_replication

echo "== host $(hostname) | class=$CLS gen=$GEN clim=$CLIM | cpus=${SLURM_CPUS_PER_TASK:-8} =="
$PY -u $PR/run_binomial.py --class "$CLS" --gen "$GEN" --climate "$CLIM" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"
echo "== ALLDONE_BINOM $CLS gen$GEN $CLIM =="
