#!/bin/bash
#SBATCH --job-name=sm_mix
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/sm_mix_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/sm_mix_%A_%a.err

# SEEDMIX_S1 cactus_em on each mixed-only cn_full variant.
# Uses the (rebuilding) cn_var_v3qc_v2 with the called mask for proper AF projection.
set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq

case $SLURM_ARRAY_TASK_ID in
  0) TAG=mixedloose   ;;
  1) TAG=mixedstrict  ;;
  2) TAG=mixedconserv ;;
  *) echo unknown; exit 1 ;;
esac

OUT_DIR=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/v3qc_v2_${TAG}_chr1
mkdir -p $OUT_DIR

READS=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix
R1=$READS/S1-1.1_P.fq.gz
R2=$READS/S1-1.2_P.fq.gz

echo "[$(date)] SEEDMIX_S1 on cn_full_v3qc_v2_${TAG} Chr1"
/usr/bin/time -v /home/tbellagio/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v3qc_v2_${TAG}/cn \
    --cn-var       data/cn_var_231_v3qc_v2.cn_var.npz \
    --cn-var-meta  data/cn_var_231_v3qc_v2.meta.npz \
    --cn-var-called data/cn_var_231_v3qc_v2.cn_var_called.npz \
    --reads $R1 $R2 \
    --sample SEEDMIX_S1 \
    --out $OUT_DIR/SEEDMIX_S1.tsv \
    --threads 8 \
    --chroms Chr1

echo "[$(date)] DONE $TAG"
ls -lh $OUT_DIR/
