#!/bin/bash
#SBATCH --job-name=kmidx_pgref
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/kmer_index/logs/pg_reference_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/kmer_index/logs/pg_reference_%j.err
set -euo pipefail

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh
conda activate pangenie

WORK=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/kmer_index

PanGenie-index \
    -r $WORK/data/test_ref_chr1.fa \
    -v $WORK/data/test_chr1_2mb.diploid.vcf \
    -o $WORK/pg_reference/pgref \
    -t 4 -k 31 \
    -e 100000000

echo "[$(date +%H:%M:%S)] DONE"
ls -lh $WORK/pg_reference/
