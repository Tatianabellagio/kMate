#!/bin/bash
#SBATCH --job-name=p231_mb
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=logs/07c_mb_%j.out
#SBATCH --error=logs/07c_mb_%j.err

# =============================================================================
# benchmarks/p231 front-runner: filt2 kmer_pa + GLOBAL EM (normalize=per_founder,
# the 2026-07-06 Kf_w fix), projected through ONE arch3 var_pa arm (atomized |
# raw). The EM h is identical across arms (same reads, same kmer_pa); only the
# projection var_pa differs.
#
# Usage: sbatch 07c_run_kmate_filt2_mb_p231.sh REGIME CNVAR [WEIGHT]
#   REGIME = n50_g0 n231_g0 n50_g1 n231_g1 n50_g3 n50_g3_dom500
#   CNVAR  = atomized | raw
#   WEIGHT = uniform (default, front-runner) | inv_mb (legacy A/B baseline)
#   GLOBAL mode drops omega=1/m_b (superseded per PIPELINE_STATE.md Sec.0,
#   2026-07-06): per_founder+uniform beats per_founder+1/m_b on AF-MAE.
#   --normalize is left at its per_sample_per_chrom.py default (per_founder).
# =============================================================================
mkdir -p logs
set -euo pipefail
REGIME=${1:?Usage: REGIME CNVAR [WEIGHT]}
CNVAR=${2:?Usage: REGIME CNVAR [WEIGHT]}
WEIGHT=${3:-uniform}
[[ "$CNVAR" == "atomized" || "$CNVAR" == "raw" ]] || { echo "ERROR: CNVAR must be atomized|raw" >&2; exit 1; }
[[ "$WEIGHT" == "inv_mb" || "$WEIGHT" == "uniform" ]] || { echo "ERROR: WEIGHT must be inv_mb|uniform" >&2; exit 1; }

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

# Front-runner kmer_pa: REBUILT-from-arch3 filt2 (single-source). Falls back to
# production v3qc_v3 filt2 if the rebuilt one is absent.
CN_KMER_PREFIX=$CTRL/data/kmer_pa_p231_filt2/kmer_pa
[ -s "${CN_KMER_PREFIX}_Chr1.kmer_pa.npz" ] || CN_KMER_PREFIX=$ROOT/data/kmer_pa_231_v3qc_v3_filt2/kmer_pa

if [ "$CNVAR" = "atomized" ]; then
    CN_VAR=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1_atomized.var_pa.npz
    CN_VAR_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1_atomized.meta.npz
else
    CN_VAR=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
    CN_VAR_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz
fi

WTAG=$([[ "$WEIGHT" == "inv_mb" ]] && echo filt2mb || echo filt2u)
OUT_DIR=$CTRL/results/kmate_global_${WTAG}_${CNVAR}/${REGIME}
mkdir -p $OUT_DIR
SAMPLE=p231_${WTAG}_${CNVAR}_${REGIME}_cov${COV}_s${SEED}
OUT_TSV=$OUT_DIR/${SAMPLE}.tsv

echo "[$(date)] kMate ${WEIGHT} GLOBAL  regime=$REGIME  var_pa=$CNVAR"
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
    --block-mode global --kmer-weight $WEIGHT

echo; echo "[$(date)] DONE -- $OUT_TSV"; ls -lh $OUT_TSV
