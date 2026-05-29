#!/bin/bash
#SBATCH --job-name=corr_em
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --array=0-12
#SBATCH --requeue
#SBATCH --output=logs/corr_%A_%a.out
#SBATCH --error=logs/corr_%A_%a.err

set -euo pipefail
ROOT=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CN=$ROOT/data/cn_full_231_v3qc_v3_filt2/cn_Chr1
OUT=$ROOT/scratch/h_fixes/correlation
mkdir -p $OUT $ROOT/src/em_fixes/logs

# 5 g0 sims (counts in scratch/g0_sweep_h_test/filt2_<SIM>.counts.npy)
# + 8 SEEDMIX reps (counts in scratch/seedmix_h_test/filt2_S{1..8}.counts.npy)
# SAMPLE  COUNTS_PATH  OUT_TAG
CFG=(
  "g0_n231_rep0_rand $ROOT/scratch/g0_sweep_h_test/filt2_g0_n231_rep0_rand.counts.npy g0_n231_rep0_rand"
  "g0_n200_rep0_rand $ROOT/scratch/g0_sweep_h_test/filt2_g0_n200_rep0_rand.counts.npy g0_n200_rep0_rand"
  "g0_n50_rep0_cact  $ROOT/scratch/g0_sweep_h_test/filt2_g0_n50_rep0_cact.counts.npy  g0_n50_rep0_cact"
  "g0_n50_rep1_bal   $ROOT/scratch/g0_sweep_h_test/filt2_g0_n50_rep1_bal.counts.npy   g0_n50_rep1_bal"
  "g0_n50_rep2_pg    $ROOT/scratch/g0_sweep_h_test/filt2_g0_n50_rep2_pg.counts.npy    g0_n50_rep2_pg"
  "S1 $ROOT/scratch/seedmix_h_test/filt2_S1.counts.npy seedmix_S1"
  "S2 $ROOT/scratch/seedmix_h_test/filt2_S2.counts.npy seedmix_S2"
  "S3 $ROOT/scratch/seedmix_h_test/filt2_S3.counts.npy seedmix_S3"
  "S4 $ROOT/scratch/seedmix_h_test/filt2_S4.counts.npy seedmix_S4"
  "S5 $ROOT/scratch/seedmix_h_test/filt2_S5.counts.npy seedmix_S5"
  "S6 $ROOT/scratch/seedmix_h_test/filt2_S6.counts.npy seedmix_S6"
  "S7 $ROOT/scratch/seedmix_h_test/filt2_S7.counts.npy seedmix_S7"
  "S8 $ROOT/scratch/seedmix_h_test/filt2_S8.counts.npy seedmix_S8"
)
read SAMPLE COUNTS TAG <<< "${CFG[$SLURM_ARRAY_TASK_ID]}"

echo "[$(date)] === correlation EM: sample=$SAMPLE tag=$TAG ==="
$PY -u $ROOT/src/em_fixes/correlation.py \
    --cn-prefix $CN \
    --counts $COUNTS \
    --sample $SAMPLE \
    --out $OUT/${TAG}.corr.npz \
    --rhos 0.0,0.5,0.9,0.99
echo "DONE $(date)"
