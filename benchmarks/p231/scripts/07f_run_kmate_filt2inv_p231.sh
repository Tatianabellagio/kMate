#!/bin/bash
#SBATCH --job-name=p231_filt2inv
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=logs/07f_filt2inv_%j.out
#SBATCH --error=logs/07f_filt2inv_%j.out

# =============================================================================
# benchmarks/p231 front-runner (2026-07-07): filt2inv kmer_pa (drops ac=1
# singletons AND ac=F invariants -- matches REAL production
# data/kmer_pa_231_arch3_filt2inv) + GLOBAL EM (normalize=per_founder, the
# 2026-07-06 Kf_w fix) + uniform kmer-weight (GLOBAL mode drops omega=1/m_b;
# PIPELINE_STATE.md Sec.0). Supersedes 07c's "filt2" (singleton-only) arm as
# the apples-to-production comparison. Only uniform is run -- inv_mb is not
# the production choice for global mode, so it's not tested here (confirmed
# 2026-07-07).
#
# Usage: sbatch 07f_run_kmate_filt2inv_p231.sh REGIME CNVAR
#   REGIME = n50_g0 n231_g0 n50_g1 n231_g1 n50_g3 n50_g1_self97 n231_g1_self97
#            n50_g3_self97 n50_g3_dom500 n50_g3_dom500nr n50_g3_dom500_self97
#   CNVAR  = atomized | raw
# =============================================================================
mkdir -p logs
set -euo pipefail
REGIME=${1:?Usage: REGIME CNVAR}
CNVAR=${2:?Usage: REGIME CNVAR}
[[ "$CNVAR" == "atomized" || "$CNVAR" == "raw" ]] || { echo "ERROR: CNVAR must be atomized|raw" >&2; exit 1; }

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
    n231_g1_self97) SUBDIR="cov10_n231_g1_s42_self97_hotspots_p231_chr1" ;;
    n50_g3_self97) SUBDIR="cov10_n50_g3_s42_self97_hotspots_p231_chr1" ;;
    n50_g3_dom500) SUBDIR="cov10_n50_g3_s42_hotspots_dom500_p231_chr1" ;;
    n50_g3_dom500nr) SUBDIR="cov10_n50_g3_s42_hotspots_dom500nr_p231_chr1" ;;
    n50_g3_dom500_self97) SUBDIR="cov10_n50_g3_s42_self97_hotspots_dom500_p231_chr1" ;;
    *) echo "ERROR: unknown REGIME '$REGIME'" >&2; exit 1 ;;
esac

WORK=$CTRL/sims/$SUBDIR
READS_DIR=$WORK/reads
[ -s "$READS_DIR/r1.fq" ] || { echo "ERROR: missing $READS_DIR/r1.fq -- run 06 first" >&2; exit 1; }

CN_KMER_PREFIX=$CTRL/data/kmer_pa_p231_filt2inv/kmer_pa
[ -s "${CN_KMER_PREFIX}_Chr1.kmer_pa.npz" ] || { echo "ERROR: missing $CN_KMER_PREFIX -- run 03c_build_kmer_pa_filt2inv_p231.sh first" >&2; exit 1; }

if [ "$CNVAR" = "atomized" ]; then
    CN_VAR=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1_atomized.var_pa.npz
    CN_VAR_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1_atomized.meta.npz
else
    CN_VAR=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
    CN_VAR_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz
fi

OUT_DIR=$CTRL/results/kmate_global_filt2invu_${CNVAR}/${REGIME}
mkdir -p $OUT_DIR
SAMPLE=p231_filt2invu_${CNVAR}_${REGIME}_cov${COV}_s${SEED}
OUT_TSV=$OUT_DIR/${SAMPLE}.tsv

echo "[$(date)] kMate filt2inv+uniform GLOBAL  regime=$REGIME  var_pa=$CNVAR"
echo "  kmer_pa: $CN_KMER_PREFIX"
echo "  var_pa:  $CN_VAR"
echo "  out:     $OUT_TSV"

$PYTHON -u $DRIVER \
    --kmer-pa-prefix $CN_KMER_PREFIX \
    --var-pa $CN_VAR \
    --var-meta $CN_VAR_META \
    --reads $READS_DIR/r1.fq $READS_DIR/r2.fq \
    --sample $SAMPLE \
    --out $OUT_TSV \
    --threads 8 --chroms Chr1 \
    --block-mode global --kmer-weight uniform

echo; echo "[$(date)] DONE -- $OUT_TSV"; ls -lh $OUT_TSV
