#!/bin/bash
#SBATCH --job-name=kmidx_pgref
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=logs/pg_reference_%j.out
#SBATCH --error=logs/pg_reference_%j.err
mkdir -p logs
set -euo pipefail

source /global/home/users/tbellg/miniforge3/etc/profile.d/conda.sh
conda activate pangenie

WORK=/global/scratch/users/tbellg/hapfire_sv/kmer_index

PanGenie-index \
    -r $WORK/data/test_ref_chr1.fa \
    -v $WORK/data/test_chr1_2mb.diploid.vcf \
    -o $WORK/pg_reference/pgref \
    -t 4 -k 31 \
    -e 100000000

echo "[$(date +%H:%M:%S)] DONE"
ls -lh $WORK/pg_reference/
