#!/bin/bash
#SBATCH --job-name=subsample_h
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --output=logs/subsample_h_%j.out
#SBATCH --error=logs/subsample_h_%j.err

set -euo pipefail
ROOT=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
OUT_DIR=$ROOT/scratch/subsample_h_test
mkdir -p $OUT_DIR $ROOT/tests/logs
mkdir -p $ROOT/data/kmer_pa_231_v3qc_v3_subsampMedian_refilt2
mkdir -p $ROOT/data/kmer_pa_231_v3qc_v3_filt2_subsampMedian

READS1=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix_trimdedup/SEEDMIX_S1_1.dedup.fq.gz
READS2=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix_trimdedup/SEEDMIX_S1_2.dedup.fq.gz

echo "=== 1. Subsample raw v3qc_v3 → median-per-stratum + refilt2 ==="
$PY -u $ROOT/src/archive/build_subsampled_cn.py \
    --in-kmer_pa $ROOT/data/kmer_pa_231_v3qc_v3/kmer_pa_Chr1.kmer_pa.npz \
    --in-meta $ROOT/data/kmer_pa_231_v3qc_v3/kmer_pa_Chr1.meta.npz \
    --out-prefix $ROOT/data/kmer_pa_231_v3qc_v3_subsampMedian_refilt2/kmer_pa_Chr1 \
    --target median --seed 42 --refilt2

echo
echo "=== 2. Subsample filt2 v3qc_v3 → median-per-stratum + refilt2 ==="
$PY -u $ROOT/src/archive/build_subsampled_cn.py \
    --in-kmer_pa $ROOT/data/kmer_pa_231_v3qc_v3_filt2/kmer_pa_Chr1.kmer_pa.npz \
    --in-meta $ROOT/data/kmer_pa_231_v3qc_v3_filt2/kmer_pa_Chr1.meta.npz \
    --out-prefix $ROOT/data/kmer_pa_231_v3qc_v3_filt2_subsampMedian/kmer_pa_Chr1 \
    --target median --seed 42 --refilt2

echo
echo "=== 3. EM h-only on subsampled-raw kmer_pa (α=0, no shape weighting) ==="
$PY -u $ROOT/src/sweep_shape_norm_h_only.py \
    --kmer_pa-prefix $ROOT/data/kmer_pa_231_v3qc_v3_subsampMedian_refilt2/kmer_pa_Chr1 \
    --reads $READS1 $READS2 \
    --sample SEEDMIX_S1 \
    --out-prefix $OUT_DIR/subsampRaw_SEEDMIX_S1_chr1 \
    --alphas 0 --threads 8 \
    --counts-cache $OUT_DIR/subsampRaw_SEEDMIX_S1_chr1.counts.npy

echo
echo "=== 4. EM h-only on subsampled-filt2 kmer_pa (α=0) ==="
$PY -u $ROOT/src/sweep_shape_norm_h_only.py \
    --kmer_pa-prefix $ROOT/data/kmer_pa_231_v3qc_v3_filt2_subsampMedian/kmer_pa_Chr1 \
    --reads $READS1 $READS2 \
    --sample SEEDMIX_S1 \
    --out-prefix $OUT_DIR/subsampFilt2_SEEDMIX_S1_chr1 \
    --alphas 0 --threads 8 \
    --counts-cache $OUT_DIR/subsampFilt2_SEEDMIX_S1_chr1.counts.npy

echo "DONE $(date)"
