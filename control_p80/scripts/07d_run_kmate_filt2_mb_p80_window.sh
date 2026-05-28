#!/bin/bash
#SBATCH --job-name=p80_mbw
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=logs/07d_mbw_%j.out
#SBATCH --error=logs/07d_mbw_%j.err

# =============================================================================
# Front-runner test, WINDOW (10 kb / star2) mode + ω=1/m_b weighting on p80.
# Mirrors 07c (global) but --block-mode window --window-bp 10000 + the standard
# anchor/smooth args. Driver now supports --kmer-weight inv_mb in window mode
# (block_em.solve_em_per_block accepts omega per-slice).
#
# Usage: sbatch 07d_run_kmate_filt2_mb_p80_window.sh REGIME [WEIGHT]
#   REGIME = n50_g0 n200_g0 n231_g0 n50_g1 n200_g1 n231_g1 n50_g3 n50_g3_dom500
#   WEIGHT = inv_mb (default) | uniform
# =============================================================================
mkdir -p logs
set -euo pipefail
REGIME=${1:?Usage: REGIME [WEIGHT]}
WEIGHT=${2:-inv_mb}
[[ "$WEIGHT" == "inv_mb" || "$WEIGHT" == "uniform" ]] || { echo "ERROR WEIGHT" >&2; exit 1; }
[[ "$WEIGHT" == "inv_mb" ]] && WTAG="filt2mbW" || WTAG="filt2uW"
[[ "$WEIGHT" == "inv_mb" ]] && ODIR="cactus_em_window_filt2_mb" || ODIR="cactus_em_window_filt2_uniform"

CTRL=/global/scratch/users/tbellg/kmate/control_p80
PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
DRIVER=/global/scratch/users/tbellg/kmate/poolfreq/src/per_sample_per_chrom.py
COV=10; SEED=42

case "$REGIME" in
    n50_g0)         N_INDIV=50;  N_GEN=0; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n80_g0)         N_INDIV=80;  N_GEN=0; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n80_g1)         N_INDIV=80;  N_GEN=1; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n200_g0)        N_INDIV=200; N_GEN=0; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n231_g0)        N_INDIV=231; N_GEN=0; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n50_g1)         N_INDIV=50;  N_GEN=1; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n200_g1)        N_INDIV=200; N_GEN=1; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n231_g1)        N_INDIV=231; N_GEN=1; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n50_g3)         N_INDIV=50;  N_GEN=3; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n50_g3_dom500)  N_INDIV=50;  N_GEN=3; SUBDIR_TAG="hotspots_dom500_p80_chr1" ;;
    *) echo "ERROR REGIME" >&2; exit 1 ;;
esac

WORK=$CTRL/sims/cov${COV}_n${N_INDIV}_g${N_GEN}_s${SEED}_${SUBDIR_TAG}
READS_DIR=$WORK/reads
for f in $READS_DIR/r1.fq $READS_DIR/r2.fq $WORK/recomb_truth.tsv.gz; do
    [ -s "$f" ] || { echo "ERROR missing $f" >&2; exit 1; }
done

CN_KMER_PREFIX=$CTRL/data/cn_full_p80_filt2/cn
CN_VAR=$CTRL/data/cn_var_p80.cn_var.npz
CN_VAR_META=$CTRL/data/cn_var_p80.meta.npz

OUT_DIR=$CTRL/results/${ODIR}/${REGIME}
mkdir -p $OUT_DIR
SAMPLE=p80_${WTAG}_${REGIME}_cov${COV}_s${SEED}
OUT_TSV=$OUT_DIR/${SAMPLE}.tsv

echo "[$(date)] kMate filt2 WINDOW(10kb) weight=$WEIGHT --regime $REGIME"
$PYTHON -u $DRIVER \
    --cn-kmer-prefix $CN_KMER_PREFIX \
    --cn-var $CN_VAR \
    --cn-var-meta $CN_VAR_META \
    --reads $READS_DIR/r1.fq $READS_DIR/r2.fq \
    --sample $SAMPLE --out $OUT_TSV \
    --threads 8 --chroms Chr1 \
    --block-mode window --window-bp 10000 \
    --global-anchor-weight 0.3 \
    --hmm-smooth-passes 5 --hmm-smooth-alpha 0.5 \
    --kmer-weight $WEIGHT

echo; echo "[$(date)] DONE -- $OUT_TSV"; ls -lh $OUT_TSV
