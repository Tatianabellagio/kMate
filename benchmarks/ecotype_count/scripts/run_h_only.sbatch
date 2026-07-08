#!/bin/bash
#SBATCH --job-name=h_only
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=benchmarks/ecotype_count/logs/honly_%x_%j.out
#SBATCH --error=benchmarks/ecotype_count/logs/honly_%x_%j.out
set -eo pipefail
# kMate --h-only (founder mixture only, no AF projection) for an ecotype-count pool.
cd /global/scratch/users/tbellg/kmate
POOL=${1:?pool}; PANEL=${2:-p80}; T=${SLURM_CPUS_PER_TASK:-8}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
case "$PANEL" in
  p80)  KP=benchmarks/p80/data/kmer_pa_p80_ours_filt2inv/kmer_pa ;;
  p231) KP=data/kmer_pa_231_arch3_filt2inv/kmer_pa ;;
esac
R1=benchmarks/$PANEL/sims/$POOL/reads/r1.fq; R2=benchmarks/$PANEL/sims/$POOL/reads/r2.fq
OUT=benchmarks/ecotype_count/results/${POOL}.tsv
echo "[$(date)] kMate --h-only on $(hostname): $POOL"
# GLOBAL mode: drop omega=1/m_b (uniform) per PIPELINE_STATE.md Sec.0 (2026-07-06).
python -m kmate.per_sample_per_chrom --kmer-pa-prefix "$KP" --h-only \
  --reads "$R1" "$R2" --sample h_${POOL} --threads "$T" --chroms Chr1 \
  --block-mode global --kmer-weight uniform --out "$OUT"
# score -> ecotype rows
python3 benchmarks/ecotype_count/scripts/score_ecotype.py \
  --h "${OUT%.tsv}.h_per_chrom.npz" \
  --weights benchmarks/$PANEL/sims/$POOL/pool_weights.tsv \
  --pool "$POOL" --panel "$PANEL" \
  --table benchmarks/ecotype_count/results/ecotype_rows/${POOL}.tsv
echo "[$(date)] DONE $POOL"
