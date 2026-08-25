#!/bin/bash
#SBATCH --job-name=wza_binom
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/wza_binom_%x_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/wza_binom_%x_%j.out
# Canonical (deg-2 Booker) WZA on a binomial-regression per-record result.
# Usage: sbatch --job-name=wzab_<class> run_wza_binomial.sh <class> [gen] [climate]
set -euo pipefail
CLS=${1:?need class}
GEN=${2:-9}
CLIM=${3:-bio1}
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication
echo "== host $(hostname) | WZA binomial $CLS gen$GEN $CLIM =="
$PY -u $PR/run_wza.py --model binomial --class "$CLS" --gen "$GEN" --climate "$CLIM"
echo "== ALLDONE_WZA_BINOM $CLS gen$GEN $CLIM =="
