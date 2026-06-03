#!/bin/bash
#SBATCH --job-name=sm_f2acw
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=1:00:00
#SBATCH --output=logs/sm_f2acw_%j.out
#SBATCH --error=logs/sm_f2acw_%j.err
mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/hapfire_sv/poolfreq
OUT_DIR=results/seedmix_231_v3_filt2_acweight_chr1
mkdir -p $OUT_DIR
/usr/bin/time -v /global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v3_filt2/cn \
    --cn-var       data/cn_var_231_v3.cn_var.npz \
    --cn-var-meta  data/cn_var_231_v3.meta.npz \
    --reads /global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz \
            /global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz \
    --sample SEEDMIX_S1 --out $OUT_DIR/SEEDMIX_S1.tsv --threads 8 --chroms Chr1 \
    --ac-weight-counts
echo "[$(date)] DONE"
