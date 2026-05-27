#!/bin/bash
#SBATCH --job-name=sm_lowmiss_test
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/sm_lowmiss_test_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/sm_lowmiss_test_%j.err

set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq

OUT_DIR=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/seedmix_v3qc_lowmiss_test
mkdir -p $OUT_DIR

READS=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix
R1=$READS/S1-1.1_P.fq.gz
R2=$READS/S1-1.2_P.fq.gz

echo "[$(date)] SEEDMIX_S1 test: cn_full_v3qc_lowmiss (F_MISSING<=0.5) + cn_var_v3qc (Chr1 only)"
/usr/bin/time -v /home/tbellagio/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v3qc_lowmiss/cn \
    --cn-var       data/cn_var_231_v3qc.cn_var.npz \
    --cn-var-meta  data/cn_var_231_v3qc.meta.npz \
    --reads $R1 $R2 \
    --sample SEEDMIX_S1 \
    --out $OUT_DIR/SEEDMIX_S1.tsv \
    --threads 8 \
    --chroms Chr1

echo "[$(date)] DONE"
ls -lh $OUT_DIR/
