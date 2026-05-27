#!/bin/bash
#SBATCH --job-name=sm_v3qc_v2_chr1
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/sm_v3qc_v2_chr1_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/sm_v3qc_v2_chr1_%j.err

# Test SEEDMIX_S1 cactus_em with the corrected v3qc_v2 cn_full (multi-allelic fix + missing-as-N).
set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq
OUT_DIR=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/seedmix_v3qc_v2_test
mkdir -p $OUT_DIR

READS=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix
R1=$READS/S1-1.1_P.fq.gz
R2=$READS/S1-1.2_P.fq.gz

CN_VAR=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_var_231_v3qc_v2.cn_var.npz
CN_VAR_META=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_var_231_v3qc_v2.meta.npz
[ -s "$CN_VAR" ] || { echo "ERROR: cn_var_v3qc_v2 not built yet"; exit 1; }

echo "[$(date)] SEEDMIX_S1 v3qc_v2 (Chr1)"
/usr/bin/time -v /home/tbellagio/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v3qc_v2/cn \
    --cn-var       $CN_VAR \
    --cn-var-meta  $CN_VAR_META \
    --reads $R1 $R2 \
    --sample SEEDMIX_S1 \
    --out $OUT_DIR/SEEDMIX_S1.tsv \
    --threads 8 \
    --chroms Chr1

echo "[$(date)] DONE"
