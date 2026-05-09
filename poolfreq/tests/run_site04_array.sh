#!/bin/bash
#SBATCH --job-name=site04_cem
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=4:00:00
#SBATCH --array=1-57%12
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/site04_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/site04_%A_%a.err

# Run cactus_em on one site04 sample per array task. 24-task concurrency.
# Reads sample N (1-indexed) from data/site04_manifest.tsv.
set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq
mkdir -p results/site04_231_v2 tests/logs

MANIFEST=data/site04_manifest.tsv
LINE=$((SLURM_ARRAY_TASK_ID + 1))   # +1 to skip header
ROW=$(awk -F'\t' -v n=$LINE 'NR==n' $MANIFEST)
SAMPLE=$(echo "$ROW" | cut -f1)
R1=$(echo "$ROW" | cut -f2)
R2=$(echo "$ROW" | cut -f3)

OUT=results/site04_231_v2/${SAMPLE}.tsv
if [ -s "$OUT" ]; then
    echo "[$(date)] $SAMPLE already done — skipping"
    exit 0
fi

echo "[$(date)] task=${SLURM_ARRAY_TASK_ID}  sample=${SAMPLE}"
echo "  R1: $R1"
echo "  R2: $R2"

/home/tbellagio/miniforge3/envs/hapfm/bin/python -u src/per_sample_driver.py \
    --cn-kmer-prefix data/cn_full_231_v2/cn \
    --cn-var data/cn_var_231_v2.cn_var.npz \
    --cn-var-meta data/cn_var_231_v2.meta.npz \
    --reads $R1 $R2 \
    --sample $SAMPLE \
    --out $OUT \
    --threads 4 \
    --block-mode window

echo "[$(date)] DONE — $OUT"
