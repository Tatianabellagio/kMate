#!/bin/bash
#SBATCH --job-name=shape_norm
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/shape_norm_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/shape_norm_%j.err

set -euo pipefail

ROOT=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
OUT_DIR=$ROOT/scratch/shape_norm_sweep_filt2
mkdir -p $OUT_DIR $ROOT/poolfreq/tests/logs

$PY $ROOT/poolfreq/src/sweep_shape_norm_h_only.py \
    --cn-prefix $ROOT/poolfreq/data/cn_full_231_v3qc_v3_filt2/cn_Chr1 \
    --reads /carnegie/nobackup/scratch/tbellagio/pang/grenenet_reads/seed_mix_trimdedup/SEEDMIX_S1_1.dedup.fq.gz \
            /carnegie/nobackup/scratch/tbellagio/pang/grenenet_reads/seed_mix_trimdedup/SEEDMIX_S1_2.dedup.fq.gz \
    --sample SEEDMIX_S1 \
    --out-prefix $OUT_DIR/SEEDMIX_S1_chr1 \
    --alphas 0,0.3,0.5,1.0,1.5,2.0,3.0 \
    --ac-threshold 4 \
    --target median \
    --threads 8 \
    --counts-cache $OUT_DIR/SEEDMIX_S1_chr1.counts.npy

echo "DONE $(date)"
