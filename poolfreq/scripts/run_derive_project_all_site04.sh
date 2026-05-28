#!/bin/bash
#SBATCH --job-name=hf_derive_proj
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=1-57
#SBATCH --requeue
#SBATCH --output=logs/hf_derive_proj_%A_%a.out
#SBATCH --error=logs/hf_derive_proj_%A_%a.err

# Per-sample: derive per-block ecotype freqs from hapFIRE outputs, then
# project through cn_var_231_v2 → per-record AFs (incl. SVs).
set -uo pipefail
mkdir -p /global/scratch/users/tbellg/kmate/poolfreq/scripts/logs
mkdir -p /global/scratch/users/tbellg/kmate/poolfreq/results/site04_hapfire_perblock_h
mkdir -p /global/scratch/users/tbellg/kmate/poolfreq/results/site04_hapfire_perblock_proj

cd /global/scratch/users/tbellg/kmate/poolfreq

LINE=$((SLURM_ARRAY_TASK_ID + 1))   # +1 to skip header
SAMPLE=$(awk -F'\t' -v n=$LINE 'NR==n {print $1}' data/site04_manifest.tsv)
[[ -z "$SAMPLE" ]] && { echo "no sample on line $LINE"; exit 1; }

H_OUT=results/site04_hapfire_perblock_h/${SAMPLE}.npz
PROJ_OUT=results/site04_hapfire_perblock_proj/${SAMPLE}.tsv

if [ -s "$PROJ_OUT" ]; then
    echo "[$(date)] $SAMPLE projection already exists — skipping"
    exit 0
fi

PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

echo "[$(date)] $SAMPLE : derive per-block h"
$PYTHON scripts/derive_hapfire_perblock_h.py \
    --block-index data/hapfire_block_index.npz \
    --uniq-haplo-tsv /global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/samples/haplotype_frequency/${SAMPLE}_unique_haplotype_frequency.txt \
    --out $H_OUT

echo "[$(date)] $SAMPLE : project to per-record AFs"
$PYTHON scripts/project_hapfire_perblock_to_records.py \
    --perblock-h $H_OUT \
    --out $PROJ_OUT

echo "[$(date)] DONE — $PROJ_OUT"
