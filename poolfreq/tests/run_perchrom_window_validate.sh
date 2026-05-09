#!/bin/bash
#SBATCH --job-name=perchrom_win
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/perchrom_win_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/perchrom_win_%j.err

# Validate per_sample_per_chrom.py in window-mode against the genome-wide
# window-mode result we already have for MLFH040120180306 (site04 task 1).
# Output goes to a separate dir so we can diff against site04_231_v2/.
set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq
mkdir -p results/site04_231_v2_perchrom

SAMPLE=MLFH040120180306
R1=/home/tbellagio/scratch/pang/grenenet_reads/grenenet-phase1/trimmed_dedup/${SAMPLE}_1.fq.gz
R2=/home/tbellagio/scratch/pang/grenenet_reads/grenenet-phase1/trimmed_dedup/${SAMPLE}_2.fq.gz

/usr/bin/time -v /home/tbellagio/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --cn-kmer-prefix data/cn_full_231_v2/cn \
    --cn-var data/cn_var_231_v2.cn_var.npz \
    --cn-var-meta data/cn_var_231_v2.meta.npz \
    --reads $R1 $R2 \
    --sample $SAMPLE \
    --out results/site04_231_v2_perchrom/${SAMPLE}.tsv \
    --threads 8 \
    --block-mode window \
    --window-bp 200000

echo "[$(date)] DONE"
ls -lh results/site04_231_v2_perchrom/
