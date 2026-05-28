#!/bin/bash
#SBATCH --job-name=sm_rn_full
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH --output=logs/sm_rn_full_%j.out
#SBATCH --error=logs/sm_rn_full_%j.err

# Full-genome (Chr1-5) SEEDMIX_S1 cactus_em on row-normalized cn_full_231_v3qc_v2.
# Confirms whether the Chr1 8.84x h-bias result holds at full-genome scale.
mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate

OUT_DIR=/global/scratch/users/tbellg/kmate/scratch/seedmix_v3qc_v2_rownorm_fullgenome
mkdir -p $OUT_DIR

READS=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix
R1=$READS/S1-1.1_P.fq.gz
R2=$READS/S1-1.2_P.fq.gz

echo "[$(date)] SEEDMIX_S1 cactus_em on row-normalized cn_full_v3qc_v2 (full genome)"
/usr/bin/time -v /global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
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
