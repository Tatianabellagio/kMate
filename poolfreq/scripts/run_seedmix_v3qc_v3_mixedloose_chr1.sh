#!/bin/bash
#SBATCH --job-name=sm_v3_ml
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/sm_v3_ml_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/sm_v3_ml_%j.err
# SEEDMIX_S1 cactus_em on v3qc_v3 mixed-loose Chr1, NEW cn_var with called mask.
#
# Reads source (UPDATED 2026-05-20):
#   /home/tbellagio/scratch/pang/grenenet_reads/seed_mix_trimdedup/  ← trim+clumpify (jobB1)
#     xwu's POOL-SEQ Trimmomatic config (ILLUMINACLIP + SLIDINGWINDOW:4:20 + LEADING:5
#     TRAILING:5 MINLEN:36) + clumpify dedup (dedupe=t dupesubs=0 optical=f) matching
#     preprocess_one.sh on the founder side. Removes ~35% PCR-duplicate k-mer count
#     inflation present in the previous trim-only FASTQs at .../seed_mix/.
#   Previous (trim-only, no dedup): /home/tbellagio/scratch/pang/grenenet_reads/seed_mix/
#     S1-1.1_P.fq.gz / S1-1.2_P.fq.gz — kept as-is on disk, not deleted.
set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq
OUT_DIR=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/v3qc_v3_mixedloose_chr1
mkdir -p $OUT_DIR
READS=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix_trimdedup
R1=$READS/SEEDMIX_S1_1.dedup.fq.gz; R2=$READS/SEEDMIX_S1_2.dedup.fq.gz

/usr/bin/time -v /home/tbellagio/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v3qc_v3_mixedloose/cn \
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
