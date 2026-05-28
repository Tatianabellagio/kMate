#!/bin/bash
#SBATCH --job-name=sm_retest
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=2:00:00
#SBATCH --output=logs/sm_retest_%A_%a.out
#SBATCH --error=logs/sm_retest_%A_%a.err

# Retest all EM variants on cn_full_v3qc_v2 Chr1 with the FIXED cn_var (called-mask
# aware projection). Each array task is a different config. After cn_var rebuild
# completes, this array re-runs SEEDMIX_S1 cactus_em through the new projection.
#
# array tasks:
#   0: baseline (no rownorm)
#   1: rownorm α=0
#   2: rownorm α=1
#   3: rownorm α=2
#   4: rownorm α=3
#   5: rownorm α=5
#   6: filt2+rownorm  (uses pre-built data/cn_full_231_v3qc_v2_filt2_rownorm/)

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate/poolfreq

PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
SCRATCH=/global/scratch/users/tbellg/kmate/scratch
READS=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix
R1=$READS/S1-1.1_P.fq.gz
R2=$READS/S1-1.2_P.fq.gz
CV=data/cn_var_231_v3qc_v2.cn_var.npz
META=data/cn_var_231_v3qc_v2.meta.npz
CALLED=data/cn_var_231_v3qc_v2.cn_var_called.npz

case $SLURM_ARRAY_TASK_ID in
  0) CN=data/cn_full_231_v3qc_v2/cn;        OUT=$SCRATCH/v3qc_v2_FIXED_baseline_chr1;        EXTRA="" ;;
  1) CN=data/cn_full_231_v3qc_v2/cn;        OUT=$SCRATCH/v3qc_v2_FIXED_kfa0p0_chr1;          EXTRA="--row-normalize-cn --kf-correction-alpha 0.0" ;;
  2) CN=data/cn_full_231_v3qc_v2/cn;        OUT=$SCRATCH/v3qc_v2_FIXED_kfa1p0_chr1;          EXTRA="--row-normalize-cn --kf-correction-alpha 1.0" ;;
  3) CN=data/cn_full_231_v3qc_v2/cn;        OUT=$SCRATCH/v3qc_v2_FIXED_kfa2p0_chr1;          EXTRA="--row-normalize-cn --kf-correction-alpha 2.0" ;;
  4) CN=data/cn_full_231_v3qc_v2/cn;        OUT=$SCRATCH/v3qc_v2_FIXED_kfa3p0_chr1;          EXTRA="--row-normalize-cn --kf-correction-alpha 3.0" ;;
  5) CN=data/cn_full_231_v3qc_v2/cn;        OUT=$SCRATCH/v3qc_v2_FIXED_kfa5p0_chr1;          EXTRA="--row-normalize-cn --kf-correction-alpha 5.0" ;;
  6) CN=data/cn_full_231_v3qc_v2_filt2_rownorm/cn; OUT=$SCRATCH/v3qc_v2_FIXED_filt2rn_chr1;  EXTRA="" ;;
  *) echo "Unknown task $SLURM_ARRAY_TASK_ID"; exit 1 ;;
esac

mkdir -p $OUT

echo "[$(date)] task=$SLURM_ARRAY_TASK_ID  CN=$CN  OUT=$OUT  EXTRA='$EXTRA'"
/usr/bin/time -v $PY -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix $CN \
    --cn-var       $CV \
    --cn-var-meta  $META \
    --cn-var-called $CALLED \
    --reads $R1 $R2 \
    --sample SEEDMIX_S1 \
    --out $OUT/SEEDMIX_S1.tsv \
    --threads 8 \
    --chroms Chr1 \
    $EXTRA

echo "[$(date)] DONE task=$SLURM_ARRAY_TASK_ID"
ls -lh $OUT/
