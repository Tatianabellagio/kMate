#!/bin/bash
#SBATCH --job-name=sm_v3_f2
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/sm_v3_f2_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/sm_v3_f2_%j.err
# SEEDMIX_S1 cactus_em on v3qc_v3 filt2 (drop ac=1 only) Chr1.
# Same cn_var + reads as the mixed-loose/strict/conserv runs — only cn_full changes.
set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq
OUT_DIR=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/v3qc_v3_filt2_chr1
mkdir -p $OUT_DIR
READS=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix_trimdedup
R1=$READS/SEEDMIX_S1_1.dedup.fq.gz; R2=$READS/SEEDMIX_S1_2.dedup.fq.gz

/usr/bin/time -v /home/tbellagio/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v3qc_v3_filt2/cn \
    --cn-var        data/cn_var_231_v3qc_v3.cn_var.npz \
    --cn-var-meta   data/cn_var_231_v3qc_v3.meta.npz \
    --cn-var-called data/cn_var_231_v3qc_v3.cn_var_called.npz \
    --reads $R1 $R2 \
    --sample SEEDMIX_S1 \
    --out $OUT_DIR/SEEDMIX_S1.tsv \
    --threads 8 \
    --chroms Chr1
echo "[$(date)] DONE"
ls -lh $OUT_DIR/
