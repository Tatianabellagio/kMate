#!/bin/bash
#SBATCH --job-name=p231_mbw
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=8:00:00
#SBATCH --output=logs/07d_mbw_%j.out
#SBATCH --error=logs/07d_mbw_%j.err

# =============================================================================
# benchmarks/p231 WINDOW-mode kMate (the recombinant-pool estimator): filt2
# kmer_pa + per-window EM + HMM smoothing + 1/m_b, projected through one arch3
# var_pa arm. Mirrors 07c (global) but --block-mode window + the star2 recipe
# (10 kb windows, anchor 0.3, 5 smooth passes). Output -> kmate_window_*.
#
# Usage: sbatch 07d_run_kmate_window_p231.sh REGIME CNVAR [WEIGHT]
#   REGIME = n50_g0 n231_g0 n50_g1 n231_g1 n50_g3 n50_g3_dom500
#            n50_g1_self97 n50_g3_self97 n231_g1_self97 n50_g3_dom500_self97
#   CNVAR  = atomized | raw      WEIGHT = inv_mb (default) | uniform
# =============================================================================
mkdir -p logs
set -euo pipefail
REGIME=${1:?Usage: REGIME CNVAR [WEIGHT]}
CNVAR=${2:?Usage: REGIME CNVAR [WEIGHT]}
WEIGHT=${3:-uniform}   # default uniform (per_founder makes inv_mb redundant); pass inv_mb to opt in
[[ "$CNVAR" == "atomized" || "$CNVAR" == "raw" ]] || { echo "ERROR: CNVAR atomized|raw" >&2; exit 1; }

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
    n50_g1_self97) SUBDIR="cov10_n50_g1_s42_self97_hotspots_p231_chr1" ;;
    n50_g3_self97) SUBDIR="cov10_n50_g3_s42_self97_hotspots_p231_chr1" ;;
    n231_g1_self97) SUBDIR="cov10_n231_g1_s42_self97_hotspots_p231_chr1" ;;
    n50_g3_dom500) SUBDIR="cov10_n50_g3_s42_hotspots_dom500_p231_chr1" ;;
    n50_g3_dom500_self97) SUBDIR="cov10_n50_g3_s42_self97_hotspots_dom500_p231_chr1" ;;
    *) echo "ERROR: unknown REGIME '$REGIME'" >&2; exit 1 ;;
esac

WORK=$CTRL/sims/$SUBDIR
READS_DIR=$WORK/reads
[ -s "$READS_DIR/r1.fq" ] || { echo "ERROR: missing $READS_DIR/r1.fq -- run 06 first" >&2; exit 1; }

CN_KMER_PREFIX=$CTRL/data/kmer_pa_p231_filt2/kmer_pa
if [ "$CNVAR" = "atomized" ]; then
    CN_VAR=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1_atomized.var_pa.npz
    CN_VAR_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1_atomized.meta.npz
else
    CN_VAR=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
    CN_VAR_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz
fi

WTAG=$([[ "$WEIGHT" == "inv_mb" ]] && echo filt2mb || echo filt2u)
OUT_DIR=$CTRL/results/kmate_window_${WTAG}_${CNVAR}/${REGIME}
mkdir -p $OUT_DIR
SAMPLE=p231_${WTAG}_${CNVAR}_${REGIME}_cov${COV}_s${SEED}
OUT_TSV=$OUT_DIR/${SAMPLE}.tsv

echo "[$(date)] kMate ${WEIGHT} WINDOW  regime=$REGIME  var_pa=$CNVAR  -> $OUT_TSV"
$PYTHON -u $DRIVER \
    --kmer-pa-prefix $CN_KMER_PREFIX \
    --var-pa $CN_VAR --var-meta $CN_VAR_META \
    --reads $READS_DIR/r1.fq $READS_DIR/r2.fq \
    --sample $SAMPLE --out $OUT_TSV \
    --threads 8 --chroms Chr1 \
    --block-mode window --window-bp 10000 \
    --global-anchor-weight 0.3 \
    --hmm-smooth-passes 5 --hmm-smooth-alpha 0.5 \
    --kmer-weight $WEIGHT
echo "[$(date)] DONE -- $OUT_TSV"; ls -lh $OUT_TSV