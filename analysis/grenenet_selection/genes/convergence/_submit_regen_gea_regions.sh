#!/bin/bash
# Re-annotate the GEA universe under the standardised rule (PROM_BP=1000, TE above promoter),
# then rebuild the hit landscape that carries those labels.
#SBATCH --job-name=regen_gea_reg
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/regen_gea_reg_%j.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "[$(date)] annotate_gea_universe"
$PY annotate_gea_universe.py
echo "[$(date)] gea_hit_landscape_data"
$PY gea_hit_landscape_data.py
echo "[$(date)] DONE"
