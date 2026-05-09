#!/bin/bash
#SBATCH --job-name=sm_perchrom_v3
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/sm_perchrom_v3_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/sm_perchrom_v3_%j.err

# v3 panel test on SEEDMIX_S1 — slope should sit at ~1.0 (no Beagle bias),
# vs v2's 1.43× scale that needed post-hoc calibration.
# Default --block-mode global — SEEDMIX is homogeneous F0 (no recomb).

set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq
mkdir -p results/seedmix_231_v3_perchrom

/usr/bin/time -v /home/tbellagio/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v3/cn \
    --cn-var data/cn_var_231_v3.cn_var.npz \
    --cn-var-meta data/cn_var_231_v3.meta.npz \
    --reads /home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz \
            /home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz \
    --sample SEEDMIX_S1 \
    --out results/seedmix_231_v3_perchrom/SEEDMIX_S1.tsv \
    --threads 8

echo "[$(date)] DONE"
ls -lh results/seedmix_231_v3_perchrom/
