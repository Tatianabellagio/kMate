#!/bin/bash
#SBATCH --job-name=sm_v3qc_f2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=2:00:00
#SBATCH --output=logs/sm_v3qc_f2_%A_%a.out
#SBATCH --error=logs/sm_v3qc_f2_%A_%a.err

# cactus_em SEEDMIX_S{N} on v3qc panel — FILT2 mode (cn_full filtered to ac_k>=2).
mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate

S=${SLURM_ARRAY_TASK_ID:-1}
OUT_DIR=/global/scratch/users/tbellg/kmate/scratch/seedmix_v3qc_filt2_dedup
mkdir -p $OUT_DIR

READS_DIR=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix
R1=$READS_DIR/S${S}-1.1_P.fq.gz
R2=$READS_DIR/S${S}-1.2_P.fq.gz
[ -s "$R1" ] || { echo "ERROR: missing $R1"; exit 1; }

echo "[$(date)] SEEDMIX_S${S} v3qc (filt2)"
/usr/bin/time -v /global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v3qc_filt2/cn \
    --cn-var       data/cn_var_231_v3qc.cn_var.npz \
    --cn-var-meta  data/cn_var_231_v3qc.meta.npz \
    --reads $R1 $R2 \
    --sample SEEDMIX_S${S} \
    --out $OUT_DIR/SEEDMIX_S${S}.tsv \
    --threads 8 \
    --chroms Chr1 Chr2 Chr3 Chr4 Chr5

echo "[$(date)] DONE"
ls -lh $OUT_DIR/SEEDMIX_S${S}.tsv
