#!/bin/bash
#SBATCH --job-name=sm_kfa
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=1:30:00
#SBATCH --requeue
#SBATCH --output=logs/sm_kfa_%A_%a.out
#SBATCH --error=logs/sm_kfa_%A_%a.err

# Sweep --kf-correction-alpha on cn_full_v3qc_v2 Chr1 (no rownorm-prebuilt).
# Loads original (unnormalized) cn_full_v3qc_v2 and applies --row-normalize-cn
# at runtime (so K_f is captured pre-norm for correction). α=1.0 reproduces
# the existing rownorm-only result; α>1 increases penalty on rich-K_f founders.
#
# IMPORTANT: We must use cn_full_v3qc_v2 (NOT _rownorm), because --row-normalize-cn
# inside the EM driver captures K_f from the original cn rows.
mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate

# Map SLURM_ARRAY_TASK_ID -> alpha value
ALPHAS=(0.0 1.0 2.0 3.0 5.0)
ALPHA=${ALPHAS[$SLURM_ARRAY_TASK_ID]}
ALPHA_TAG=${ALPHA//./p}

OUT_DIR=/global/scratch/users/tbellg/kmate/scratch/seedmix_v3qc_v2_kfalpha${ALPHA_TAG}_chr1
mkdir -p $OUT_DIR

READS=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix
R1=$READS/S1-1.1_P.fq.gz
R2=$READS/S1-1.2_P.fq.gz

echo "[$(date)] SEEDMIX_S1 cactus_em on v3qc_v2 Chr1 --row-normalize-cn --kf-correction-alpha=$ALPHA"
/usr/bin/time -v /global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v3qc_v2/cn \
    --cn-var       data/cn_var_231_v3qc_v2.cn_var.npz \
    --cn-var-meta  data/cn_var_231_v3qc_v2.meta.npz \
    --reads $R1 $R2 \
    --sample SEEDMIX_S1 \
    --out $OUT_DIR/SEEDMIX_S1.tsv \
    --threads 8 \
    --chroms Chr1 \
    --row-normalize-cn \
    --kf-correction-alpha $ALPHA

echo "[$(date)] DONE α=$ALPHA"
ls -lh $OUT_DIR/
