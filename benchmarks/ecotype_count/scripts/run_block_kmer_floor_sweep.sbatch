#!/bin/bash
#SBATCH --job-name=kmer_floor_sw
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=1:30:00
#SBATCH --array=0-15
#SBATCH --output=benchmarks/ecotype_count/logs/%x_%A_%a.out
#SBATCH --error=benchmarks/ecotype_count/logs/%x_%A_%a.out
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate

# Window-size sweep to EXTEND the k-mers-per-block axis up to the plateau, on the
# few-ecotype pools where "differentiate ecotypes" is a real question. Combined
# with the 5 kb run, blocks then span ~1 .. ~50k k-mers -> the accuracy curve
# reaches the global solution and the floor becomes visible.
POOLS=(cov10_n3_g0_s42_hotspots_p80_chr1
       cov10_n5_g0_s42_hotspots_p80_chr1
       cov10_n10_g0_s42_hotspots_p80_chr1
       cov10_n50_g0_s42_hotspots_p80_chr1)
WINS=(20000 80000 320000 1280000)

P=$(( SLURM_ARRAY_TASK_ID / 4 ))
W=$(( SLURM_ARRAY_TASK_ID % 4 ))
POOL=${POOLS[$P]}; WIN=${WINS[$W]}
OUT=benchmarks/ecotype_count/results/kmer_floor/${POOL}__w${WIN}.tsv
echo "[$(date)] host=$(hostname) pool=$POOL window=$WIN"

python benchmarks/ecotype_count/scripts/bench_block_kmer_floor.py \
  --pool "$POOL" --panel p80 --chrom Chr1 \
  --window-bp "$WIN" --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --kmer-weight inv_mb --out "$OUT"
echo "[$(date)] done $POOL w$WIN"
