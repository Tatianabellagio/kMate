#!/bin/bash
#SBATCH --job-name=cc_filter
#SBATCH --partition=bse
#SBATCH --cpus-per-task=2
#SBATCH --mem=192G
#SBATCH --time=0:30:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/block_haplotype_cn/cc_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/block_haplotype_cn/cc_%j.err

set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv

PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python
INPUT=${1:?input npz}
MIN_CC=${MIN_CC:-5}
TAG="cc${MIN_CC}"
OUTPUT=$(echo $INPUT | sed "s/.npz$/_${TAG}.npz/")

echo "[$(date)] filtering $INPUT (min_cc=$MIN_CC) → $OUTPUT"
$PYTHON -u poolfreq/scripts/filter_block_haplotype_cn_min_carrier.py \
    --in $INPUT --out $OUTPUT --min-cc $MIN_CC

echo "[$(date)] DONE"
ls -lh $OUTPUT
