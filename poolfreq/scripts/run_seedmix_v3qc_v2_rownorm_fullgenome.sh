#!/bin/bash
#SBATCH --job-name=sm_rn_full
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/sm_rn_full_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/sm_rn_full_%j.err

# Full-genome (Chr1-5) SEEDMIX_S1 cactus_em on row-normalized cn_full_231_v3qc_v2.
# Confirms whether the Chr1 8.84x h-bias result holds at full-genome scale.
set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq

OUT_DIR=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/seedmix_v3qc_v2_rownorm_fullgenome
mkdir -p $OUT_DIR

READS=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix
R1=$READS/S1-1.1_P.fq.gz
R2=$READS/S1-1.2_P.fq.gz

echo "[$(date)] SEEDMIX_S1 cactus_em on row-normalized cn_full_v3qc_v2 (full genome)"
/usr/bin/time -v /home/tbellagio/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v3qc_v2_rownorm/cn \
    --cn-var       data/cn_var_231_v3qc_v2.cn_var.npz \
    --cn-var-meta  data/cn_var_231_v3qc_v2.meta.npz \
    --reads $R1 $R2 \
    --sample SEEDMIX_S1 \
    --out $OUT_DIR/SEEDMIX_S1.tsv \
    --threads 8 \
    --chroms Chr1 Chr2 Chr3 Chr4 Chr5

echo "[$(date)] DONE"
ls -lh $OUT_DIR/
