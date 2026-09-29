#!/bin/bash
#SBATCH --job-name=p231_himb
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --output=logs/09_himb_%j.out
#SBATCH --error=logs/09_himb_%j.out

# =============================================================================
# h-imbalance grid, PRODUCTION --unit chrom (per_founder normalize, default).
# Two-row version: raw (no filter) vs filt2inv (production filter) -- drops the
# old third row (filt2 + bubble weighting omega=1/m_b), since per-founder
# M-step normalization supersedes that weighting (2026-07-06 fix). h-only, no
# var-pa needed (matches notebooks/_build_h_imbalance_nb.py's chrom_*.h_per_chrom.npz
# consumer contract).
#
# Usage: sbatch 09_run_h_imbalance_chrom_p231.sh REGIME ARM
#   REGIME = n231_g0 | n50_cactheavy | n50_pgheavy
#   ARM    = raw | filt2inv
# =============================================================================
mkdir -p logs
set -euo pipefail
REGIME=${1:?Usage: REGIME ARM}
ARM=${2:?Usage: REGIME ARM}
ROOT=/global/scratch/users/tbellg/kmate
CTRL=$ROOT/benchmarks/p231
PYTHON=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
DRIVER=$ROOT/src/per_sample_per_chrom.py
COV=10; SEED=42

case "$REGIME" in
    n231_g0)       SUBDIR="cov10_n231_g0_s42_hotspots_p231_chr1" ;;
    n50_cactheavy) SUBDIR="cov10_n50_cactheavy_s42_hotspots_p231_chr1" ;;
    n50_pgheavy)   SUBDIR="cov10_n50_pgheavy_s42_hotspots_p231_chr1" ;;
    *) echo "ERROR: unknown REGIME '$REGIME'" >&2; exit 1 ;;
esac
case "$ARM" in
    raw)      KMER=$CTRL/data/kmer_pa_p231/kmer_pa ;;
    filt2inv) KMER=$CTRL/data/kmer_pa_p231_filt2inv/kmer_pa ;;
    *) echo "ERROR: ARM must be raw|filt2inv" >&2; exit 1 ;;
esac

READS_DIR=$CTRL/sims/$SUBDIR/reads
[ -s "$READS_DIR/r1.fq" ] || { echo "ERROR: missing $READS_DIR/r1.fq" >&2; exit 1; }
[ -s "${KMER}_Chr1.kmer_pa.npz" ] || { echo "ERROR: missing $KMER" >&2; exit 1; }

OUT_DIR=$CTRL/results/kmate_chrom_himbalance_${ARM}/${REGIME}
mkdir -p $OUT_DIR
SAMPLE=p231_chrom_himbalance_${ARM}_${REGIME}_cov${COV}_s${SEED}
OUT_TSV=$OUT_DIR/${SAMPLE}.tsv

echo "[$(date)] kMate --unit chrom  h-only  regime=$REGIME  arm=$ARM"
echo "  kmer_pa: $KMER"

$PYTHON -u $DRIVER \
    --kmer-pa-prefix $KMER \
    --h-only \
    --reads $READS_DIR/r1.fq $READS_DIR/r2.fq \
    --sample $SAMPLE \
    --out $OUT_TSV \
    --threads 4 --chroms Chr1 \
    --unit chrom --kmer-weight uniform

echo; echo "[$(date)] DONE — ${OUT_TSV%.tsv}.h_per_chrom.npz"
