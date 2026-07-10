#!/bin/bash
#SBATCH --job-name=block_nofb_af
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:30:00
#SBATCH --array=0-1
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/%x_%A_%a.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/%x_%A_%a.out
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
export PYTHONPATH=/global/scratch/users/tbellg/kmate/src:$PYTHONPATH

# Block-mode AF projection with the GLOBAL FALLBACK DISABLED (min-kmers-per-block 1).
# task 0: pure local  -> anchor 0, no smoothing  (zero global influence anywhere)
# task 1: no-fallback but anchored/smoothed (production niceties, still no fallback)
POOL=cov10_n80_g0_s42_hotspots_p80_chr1
D=benchmarks/p80
WORK=$D/sims/$POOL
KP=$D/data/kmer_pa_p80_filt2/kmer_pa
VAR=$D/data/var_pa_p80.var_pa.npz
META=$D/data/var_pa_p80.meta.npz
CALLED=$D/data/var_pa_p80.var_called.npz
OUTDIR=benchmarks/accuracy_vs_competitors/work
mkdir -p "$OUTDIR"

if [ "$SLURM_ARRAY_TASK_ID" = "0" ]; then
  TAG=blocknofb_purelocal; ANCHOR=0.0; SMOOTH=0
else
  TAG=blocknofb_anchored;  ANCHOR=0.3; SMOOTH=5
fi
OUT=$OUTDIR/${POOL}__${TAG}.tsv
echo "[$(date)] host=$(hostname) $TAG anchor=$ANCHOR smooth=$SMOOTH window=10kb min-kmers=1 (no fallback)"

python -m kmate.per_sample_per_chrom \
  --kmer-pa-prefix "$KP" --var-pa "$VAR" --var-meta "$META" --var-called "$CALLED" \
  --reads "$WORK/reads/r1.fq" "$WORK/reads/r2.fq" \
  --sample "$TAG" --out "$OUT" \
  --threads "${SLURM_CPUS_PER_TASK:-8}" --chroms Chr1 \
  --block-mode window --window-bp 10000 --kmer-weight inv_mb \
  --min-kmers-per-block 1 \
  --global-anchor-weight "$ANCHOR" --hmm-smooth-passes "$SMOOTH"
echo "[$(date)] done -> $OUT"
