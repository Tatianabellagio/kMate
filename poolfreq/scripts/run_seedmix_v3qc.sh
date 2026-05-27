#!/bin/bash
#SBATCH --job-name=sm_v3qc
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/sm_v3qc_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/sm_v3qc_%A_%a.err

# cactus_em SEEDMIX_S{N} on v3qc panel — REGULAR mode (no filt2).
# Array task ID = SEEDMIX replicate (1..8).
# Output: scratch/seedmix_v3qc_dedup/SEEDMIX_S{N}.tsv
set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq

S=${SLURM_ARRAY_TASK_ID:-1}
OUT_DIR=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/seedmix_v3qc_dedup
mkdir -p $OUT_DIR

READS_DIR=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix
R1=$READS_DIR/S${S}-1.1_P.fq.gz
R2=$READS_DIR/S${S}-1.2_P.fq.gz
[ -s "$R1" ] || { echo "ERROR: missing $R1"; exit 1; }

echo "[$(date)] SEEDMIX_S${S} v3qc (regular)"
/usr/bin/time -v /home/tbellagio/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v3qc/cn \
    --cn-var       data/cn_var_231_v3qc.cn_var.npz \
    --cn-var-meta  data/cn_var_231_v3qc.meta.npz \
    --reads $R1 $R2 \
    --sample SEEDMIX_S${S} \
    --out $OUT_DIR/SEEDMIX_S${S}.tsv \
    --threads 8 \
    --chroms Chr1 Chr2 Chr3 Chr4 Chr5

echo "[$(date)] DONE"
ls -lh $OUT_DIR/SEEDMIX_S${S}.tsv
