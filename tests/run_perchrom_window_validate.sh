#!/bin/bash
#SBATCH --job-name=perchrom_win
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --output=logs/perchrom_win_%j.out
#SBATCH --error=logs/perchrom_win_%j.err

# Validate per_sample_per_chrom.py in window-mode against the genome-wide
# window-mode result we already have for MLFH040120180306 (site04 task 1).
# Output goes to a separate dir so we can diff against site04_231_v2/.
mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p results/site04_231_v2_perchrom

SAMPLE=MLFH040120180306
R1=/global/scratch/users/tbellg/pang/grenenet_reads/grenenet-phase1/trimmed_dedup/${SAMPLE}_1.fq.gz
R2=/global/scratch/users/tbellg/pang/grenenet_reads/grenenet-phase1/trimmed_dedup/${SAMPLE}_2.fq.gz

/usr/bin/time -v /global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u src/per_sample_per_chrom.py \
    --kmer-pa-prefix data/kmer_pa_231_v2/kmer_pa \
    --var-pa data/var_pa_231_v2.var_pa.npz \
    --var-meta data/var_pa_231_v2.meta.npz \
    --reads $R1 $R2 \
    --sample $SAMPLE \
    --out results/site04_231_v2_perchrom/${SAMPLE}.tsv \
    --threads 8 \
    --block-mode window \
    --window-bp 200000

echo "[$(date)] DONE"
ls -lh results/site04_231_v2_perchrom/
