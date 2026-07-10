#!/bin/bash
#SBATCH --job-name=cov_sim_km
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/covsim_%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/covsim_%x_%j.err
set -eo pipefail
# Generate a p80 sim pool at arbitrary coverage + run kMate on it.
# Usage: sbatch run_cov_sim_kmate.sbatch <REGIME> <COV> <N_INDIV> <N_GEN>
#   e.g. ... n80_g0 50 80 0
cd /global/scratch/users/tbellg/kmate/benchmarks/p80/scripts
REG=${1:?regime e.g. n80_g0}; export COV=${2:?coverage e.g. 50}
N=${3:?N_INDIV}; G=${4:?N_GEN}
echo "[$(date)] cov${COV} sim+kMate for $REG (N=$N G=$G) on $(hostname)"
echo "[$(date)] === STAGE: simulate reads at ${COV}x ==="
bash 06_run_sim_p80_covarg.sh "$N" "$G"
echo "[$(date)] === STAGE: kMate (filt2 inv_mb, global) ==="
bash 07c_run_kmate_filt2_mb_p80_covarg.sh "$REG" inv_mb
echo "[$(date)] DONE cov${COV} $REG"