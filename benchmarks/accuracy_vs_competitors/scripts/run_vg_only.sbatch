#!/bin/bash
#SBATCH --job-name=vg_only
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/vgonly_%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/vgonly_%x_%j.out
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
POOL=${1:?pool}; PANEL=${2:-p80}; PRESET=${3:-}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vg
/usr/bin/time -v bash benchmarks/accuracy_vs_competitors/scripts/run_vg.sh "$POOL" "$PANEL" "${SLURM_CPUS_PER_TASK:-8}" "$PRESET"
