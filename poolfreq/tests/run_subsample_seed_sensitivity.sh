#!/bin/bash
#SBATCH --job-name=subsamp_seeds
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/subsamp_seeds_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/subsamp_seeds_%j.err

set -euo pipefail
ROOT=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
OUT_DIR=$ROOT/scratch/subsample_h_test
READS1=/carnegie/nobackup/scratch/tbellagio/pang/grenenet_reads/seed_mix_trimdedup/SEEDMIX_S1_1.dedup.fq.gz
READS2=/carnegie/nobackup/scratch/tbellagio/pang/grenenet_reads/seed_mix_trimdedup/SEEDMIX_S1_2.dedup.fq.gz

for SEED in 1 7 100; do
    echo "=== Seed $SEED: subsample raw v3qc_v3 → refilt2 → EM ==="
    CN_DIR=$ROOT/poolfreq/data/cn_full_231_v3qc_v3_subsampMedian_refilt2_seed${SEED}
    mkdir -p $CN_DIR

    $PY -u $ROOT/poolfreq/src/build_subsampled_cn.py \
        --in-cn $ROOT/poolfreq/data/cn_full_231_v3qc_v3/cn_Chr1.cn.npz \
        --in-meta $ROOT/poolfreq/data/cn_full_231_v3qc_v3/cn_Chr1.meta.npz \
        --out-prefix $CN_DIR/cn_Chr1 \
        --target median --seed $SEED --refilt2

    $PY -u $ROOT/poolfreq/src/sweep_shape_norm_h_only.py \
        --cn-prefix $CN_DIR/cn_Chr1 \
        --reads $READS1 $READS2 \
        --sample SEEDMIX_S1 \
        --out-prefix $OUT_DIR/subsampRaw_seed${SEED}_SEEDMIX_S1_chr1 \
        --alphas 0 --threads 8 \
        --counts-cache $OUT_DIR/subsampRaw_seed${SEED}_SEEDMIX_S1_chr1.counts.npy
done

echo "DONE $(date)"
