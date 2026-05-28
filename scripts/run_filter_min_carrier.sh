#!/bin/bash
#SBATCH --job-name=cc_filter
#SBATCH --account=fc_moilab
#SBATCH --partition=savio3_bigmem
#SBATCH --qos=savio_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=192G
#SBATCH --time=0:30:00
#SBATCH --output=logs/cc_%j.out
#SBATCH --error=logs/cc_%j.err

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate

PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
INPUT=${1:?input npz}
MIN_CC=${MIN_CC:-5}
TAG="cc${MIN_CC}"
OUTPUT=$(echo $INPUT | sed "s/.npz$/_${TAG}.npz/")

echo "[$(date)] filtering $INPUT (min_cc=$MIN_CC) → $OUTPUT"
$PYTHON -u scripts/filter_block_haplotype_cn_min_carrier.py \
    --in $INPUT --out $OUTPUT --min-cc $MIN_CC

echo "[$(date)] DONE"
ls -lh $OUTPUT
