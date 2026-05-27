#!/bin/bash
#SBATCH --job-name=sm_filt2
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=1:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/jf_chr1/logs/sm_filt2_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/jf_chr1/logs/sm_filt2_%j.err

set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq

OUT_DIR=results/seedmix_231_v3_filt2_chr1
mkdir -p $OUT_DIR

/usr/bin/time -v /home/tbellagio/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v3_filt2/cn \
    --cn-var       data/cn_var_231_v3.cn_var.npz \
    --cn-var-meta  data/cn_var_231_v3.meta.npz \
    --reads /home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz \
            /home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz \
    --sample SEEDMIX_S1 \
    --out $OUT_DIR/SEEDMIX_S1.tsv \
    --threads 8 \
    --chroms Chr1

echo "[$(date)] DONE"
ls -lh $OUT_DIR/
