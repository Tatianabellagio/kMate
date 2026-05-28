#!/bin/bash
#SBATCH --job-name=genuniq_filter
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=logs/genuniq_%j.out
#SBATCH --error=logs/genuniq_%j.err

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/hapfire_sv

PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
INPUT=${1:-poolfreq/data/block_haplotype_cn/chr1_full.npz}
# Output: <input_stem>_genuniq.npz
OUTPUT=$(echo $INPUT | sed 's/.npz$/_genuniq.npz/')

echo "[$(date)] filtering $INPUT → $OUTPUT"
$PYTHON -u poolfreq/scripts/filter_block_haplotype_cn_genomic_uniqueness.py \
    --in $INPUT --out $OUTPUT \
    --fastas-dir sims/visor_freqk/founder_fastas_231 \
    --chrom Chr1 --threads 2

echo "[$(date)] DONE"
ls -lh $OUTPUT
