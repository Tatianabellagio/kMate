#!/bin/bash
#SBATCH --job-name=shape_norm
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=logs/shape_norm_%j.out
#SBATCH --error=logs/shape_norm_%j.err

set -euo pipefail

ROOT=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
OUT_DIR=$ROOT/scratch/shape_norm_sweep_filt2
mkdir -p $OUT_DIR $ROOT/tests/logs

$PY $ROOT/src/sweep_shape_norm_h_only.py \
    --cn-prefix $ROOT/data/cn_full_231_v3qc_v3_filt2/cn_Chr1 \
    --reads /global/scratch/users/tbellg/pang/grenenet_reads/seed_mix_trimdedup/SEEDMIX_S1_1.dedup.fq.gz \
            /global/scratch/users/tbellg/pang/grenenet_reads/seed_mix_trimdedup/SEEDMIX_S1_2.dedup.fq.gz \
    --sample SEEDMIX_S1 \
    --out-prefix $OUT_DIR/SEEDMIX_S1_chr1 \
    --alphas 0,0.3,0.5,1.0,1.5,2.0,3.0 \
    --ac-threshold 4 \
    --target median \
    --threads 8 \
    --counts-cache $OUT_DIR/SEEDMIX_S1_chr1.counts.npy

echo "DONE $(date)"
