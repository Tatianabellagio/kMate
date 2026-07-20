#!/bin/bash
#SBATCH --job-name=cluster_mcf90
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/cluster_mcf90_%A_%a.out
#SBATCH --error=logs/cluster_mcf90_%A_%a.out
# Downstream re-cluster on the NEW min-called-frac 0.9 blocks (blocks_mcf90).
# Per chrom: block_unit_frontier -> per-block (n_eff, block_ve) + per-cluster (cluster_ve)
# on the gen9 pools. CLQ0.9 + gate2. Env vars align the founder matrix filter (MINCF=0.9)
# and block/output dir to the rebuilt map.

set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
export MINCF=0.9
export BLOCKS_DIR=analysis/grenenet_gea/blocks_mcf90

CHRS=(Chr1 Chr2 Chr3 Chr4 Chr5)
CHR=${CHRS[$SLURM_ARRAY_TASK_ID]}
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python

echo "[$(date)] $CHR : block_unit_frontier (CLQ0.9, gate2) on blocks_mcf90 (MINCF=$MINCF)"
$PY -u analysis/grenenet_gea/block_unit_frontier.py --clqcut 0.9 --chrom "$CHR" --cluster-gate 2.0
echo "[$(date)] $CHR done"
