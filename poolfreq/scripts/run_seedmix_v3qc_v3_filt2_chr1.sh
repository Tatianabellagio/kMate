#!/bin/bash
#SBATCH --job-name=sm_v3_f2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=2:00:00
#SBATCH --output=logs/sm_v3_f2_%j.out
#SBATCH --error=logs/sm_v3_f2_%j.err
# SEEDMIX_S1 cactus_em on v3qc_v3 filt2 (drop ac=1 only) Chr1.
# Same cn_var + reads as the mixed-loose/strict/conserv runs — only cn_full changes.
mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/hapfire_sv/poolfreq
OUT_DIR=/global/scratch/users/tbellg/hapfire_sv/scratch/v3qc_v3_filt2_chr1
mkdir -p $OUT_DIR
READS=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix_trimdedup
R1=$READS/SEEDMIX_S1_1.dedup.fq.gz; R2=$READS/SEEDMIX_S1_2.dedup.fq.gz

/usr/bin/time -v /global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
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
