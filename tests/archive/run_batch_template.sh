#!/bin/bash
#SBATCH --job-name=kmate_batch
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=0-49%10            # 50 array tasks, max 10 concurrent (adjust)
#SBATCH --requeue
#SBATCH --output=logs/batch_%A_%a.out
#SBATCH --error=logs/batch_%A_%a.err

# SLURM array template for processing N samples in parallel.
# Each array task processes a chunk of samples from the manifest.
#
# Configuration (env-overridable):
#   MANIFEST         — TSV with sample_id, reads_path, [reads_path2]
#   OUT_DIR          — output directory
#   CHUNK_SIZE       — samples per array task
#   CN_KMER_PREFIX   — path prefix to per-chrom kmer_pa matrices
#   CN_VAR           — path to var_pa matrix
#   CN_VAR_META      — path to var_pa meta
#
# To use:
#   1. Build manifest from your sample sources
#   2. Set --array=0-N where N = ceil(n_samples / CHUNK_SIZE) - 1
#   3. sbatch (with env overrides as needed)
#
# Examples:
#   # 82-founder run on SEEDMIX:
#   MANIFEST=data/seedmix_manifest.tsv OUT_DIR=results/seedmix_82 sbatch --array=0-0 ...
#
#   # 231-founder run on the 2,415 evolved samples (after imputation lands):
#   MANIFEST=data/sample_manifest.tsv \
#   OUT_DIR=results/per_sample_231 \
#   CN_KMER_PREFIX=data/kmer_pa_231/kmer_pa \
#   CN_VAR=data/var_pa_231.var_pa.npz \
#   CN_VAR_META=data/var_pa_231.meta.npz \
#     sbatch --array=0-48%10 tests/run_batch_template.sh

set -uo pipefail
mkdir -p /global/scratch/users/tbellg/kmate/logs
cd /global/scratch/users/tbellg/kmate

MANIFEST="${MANIFEST:-data/sample_manifest.tsv}"
OUT_DIR="${OUT_DIR:-results/per_sample}"
CHUNK_SIZE="${CHUNK_SIZE:-50}"
CN_KMER_PREFIX="${CN_KMER_PREFIX:-data/kmer_pa}"
CN_VAR="${CN_VAR:-data/var_pa_82.var_pa.npz}"
CN_VAR_META="${CN_VAR_META:-data/var_pa_82.meta.npz}"

# Extract this task's chunk
START=$((SLURM_ARRAY_TASK_ID * CHUNK_SIZE + 2))   # +2 to skip header (row 1) and 1-index
END=$((START + CHUNK_SIZE - 1))

echo "Array task $SLURM_ARRAY_TASK_ID processing manifest rows $START..$END"

# Build chunk manifest (header + rows in range)
CHUNK_MANIFEST="${OUT_DIR}/_chunks/manifest_${SLURM_ARRAY_TASK_ID}.tsv"
mkdir -p "$(dirname $CHUNK_MANIFEST)"
head -1 $MANIFEST > $CHUNK_MANIFEST
sed -n "${START},${END}p" $MANIFEST >> $CHUNK_MANIFEST

/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python src/batch_runner.py \
    --manifest $CHUNK_MANIFEST \
    --out-dir $OUT_DIR \
    --kmer-pa-prefix $CN_KMER_PREFIX \
    --var-pa $CN_VAR \
    --var-meta $CN_VAR_META \
    --threads 8 \
    --workers 1
