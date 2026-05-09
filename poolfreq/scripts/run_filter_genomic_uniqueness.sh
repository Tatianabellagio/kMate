#!/bin/bash
#SBATCH --job-name=genuniq_filter
#SBATCH --partition=bse
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/block_haplotype_cn/genuniq_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/block_haplotype_cn/genuniq_%j.err

set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv

PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python
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
