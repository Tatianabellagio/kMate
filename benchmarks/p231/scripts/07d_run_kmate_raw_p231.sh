#!/bin/bash
#SBATCH --job-name=p231_raw
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/benchmarks/p231/logs/07d_raw_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/benchmarks/p231/logs/07d_raw_%j.err
# 07d -- RAW (unfiltered kmer_pa, uniform weight) arm for the h-imbalance grid.
# Mirrors 07c but points --kmer-pa-prefix at the UNFILTERED kmer_pa_p231 and uses
# --kmer-weight uniform. Produces the "raw (no filter, no correction)" row.
# Reuses each regime's already-simulated reads (no re-sim). kmate env.
#
# Usage: bash 07d_run_kmate_raw_p231.sh REGIME   (REGIME = n231_g0 n50_g0 n50_g1 ...)
set -euo pipefail
REGIME=${1:?Usage: REGIME}
ROOT=/global/scratch/users/tbellg/kmate
CTRL=$ROOT/benchmarks/p231
PYTHON=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
DRIVER=$ROOT/src/per_sample_per_chrom.py
COV=10; SEED=42

case "$REGIME" in
    n50_g0)        SUBDIR="cov10_n50_g0_s42_hotspots_p231_chr1" ;;
    n231_g0)       SUBDIR="cov10_n231_g0_s42_hotspots_p231_chr1" ;;
    n50_g1)        SUBDIR="cov10_n50_g1_s42_hotspots_p231_chr1" ;;
    n231_g1)       SUBDIR="cov10_n231_g1_s42_hotspots_p231_chr1" ;;
    n50_g3)        SUBDIR="cov10_n50_g3_s42_hotspots_p231_chr1" ;;
    n50_g3_dom500) SUBDIR="cov10_n50_g3_s42_hotspots_dom500_p231_chr1" ;;
    *) echo "ERROR: unknown REGIME '$REGIME'" >&2; exit 1 ;;
esac

WORK=$CTRL/sims/$SUBDIR
READS_DIR=$WORK/reads
[ -s "$READS_DIR/r1.fq" ] || { echo "ERROR: missing $READS_DIR/r1.fq" >&2; exit 1; }

# RAW unfiltered kmer_pa (the "no filter" arm).
CN_KMER_PREFIX=$CTRL/data/kmer_pa_p231/kmer_pa
CN_VAR=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
CN_VAR_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz

OUT_DIR=$CTRL/results/kmate_global_raw_raw/${REGIME}
mkdir -p $OUT_DIR
SAMPLE=p231_raw_raw_${REGIME}_cov${COV}_s${SEED}
OUT_TSV=$OUT_DIR/${SAMPLE}.tsv

echo "[$(date)] kMate RAW (uniform) GLOBAL  regime=$REGIME"
echo "  kmer_pa: $CN_KMER_PREFIX  (unfiltered)"
echo "  out:     $OUT_TSV"

$PYTHON -u $DRIVER \
    --kmer-pa-prefix $CN_KMER_PREFIX \
    --var-pa $CN_VAR \
    --var-meta $CN_VAR_META \
    --reads $READS_DIR/r1.fq $READS_DIR/r2.fq \
    --sample $SAMPLE \
    --out $OUT_TSV \
    --threads 4 --chroms Chr1 \
    --block-mode global --kmer-weight uniform

echo "[$(date)] DONE -- $OUT_TSV  (h: ${OUT_TSV%.tsv}.h_per_chrom.npz)"
