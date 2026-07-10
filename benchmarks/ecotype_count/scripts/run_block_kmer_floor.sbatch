#!/bin/bash
#SBATCH --job-name=kmer_floor
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=1:30:00
#SBATCH --array=0-5
#SBATCH --output=benchmarks/ecotype_count/logs/%x_%A_%a.out
#SBATCH --error=benchmarks/ecotype_count/logs/%x_%A_%a.out
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate  # needs jellyfish (basic lacks it)

# g0 pools: mixture sweep at cov10 + a cov50 contrast for n80. (g0 => every
# block's true mixture == the global pool_weights, so per-block accuracy vs
# k-mer count is a clean floor curve.)
POOLS=(
  cov10_n3_g0_s42_hotspots_p80_chr1
  cov10_n5_g0_s42_hotspots_p80_chr1
  cov10_n10_g0_s42_hotspots_p80_chr1
  cov10_n50_g0_s42_hotspots_p80_chr1
  cov10_n80_g0_s42_hotspots_p80_chr1
  cov50_n80_g0_s42_hotspots_p80_chr1
)
POOL=${POOLS[$SLURM_ARRAY_TASK_ID]}
OUT=benchmarks/ecotype_count/results/kmer_floor/${POOL}.tsv
echo "[$(date)] host=$(hostname) task=$SLURM_ARRAY_TASK_ID pool=$POOL"

python benchmarks/ecotype_count/scripts/bench_block_kmer_floor.py \
  --pool "$POOL" --panel p80 --chrom Chr1 \
  --window-bp 5000 --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --kmer-weight inv_mb --out "$OUT"
echo "[$(date)] done $POOL"
