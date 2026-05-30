#!/bin/bash
#SBATCH --job-name=p80_mb
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/07c_mb_%j.out
#SBATCH --error=logs/07c_mb_%j.err

# =============================================================================
# Front-runner test: filt2 kmer_pa + GLOBAL EM + ω_k=1/m_b weighting.
# Mirrors 07b (filt2) but adds --kmer-weight inv_mb. GLOBAL mode only.
#   --kmer-pa-prefix .../kmer_pa_p80_filt2/kmer_pa
#   --kmer-weight inv_mb
#   outputs to results/cactus_em_global_filt2_mb/<REGIME>/p80_filt2mb_*.tsv
#
# Usage: sbatch 07c_run_kmate_filt2_mb_p80.sh REGIME [WEIGHT]
#   REGIME = n50_g0 n200_g0 n231_g0 n50_g1 n200_g1 n231_g1 n50_g3 n50_g3_dom500
#   WEIGHT = inv_mb (default, front-runner) | uniform (filt2 baseline for A/B)
# =============================================================================
mkdir -p logs
set -euo pipefail

REGIME=${1:?Usage: sbatch 07c_run_kmate_filt2_mb_p80.sh REGIME [WEIGHT]}
WEIGHT=${2:-inv_mb}
if [[ "$WEIGHT" != "inv_mb" && "$WEIGHT" != "uniform" ]]; then
    echo "ERROR: WEIGHT must be 'inv_mb' or 'uniform'; got '$WEIGHT'" >&2; exit 1
fi
[[ "$WEIGHT" == "inv_mb" ]] && WTAG="filt2mb" || WTAG="filt2u"
[[ "$WEIGHT" == "inv_mb" ]] && ODIR="cactus_em_global_filt2_mb" || ODIR="cactus_em_global_filt2_uniform"

CTRL=/global/scratch/users/tbellg/kmate/benchmarks/p80
PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
DRIVER=/global/scratch/users/tbellg/kmate/src/per_sample_per_chrom.py

COV=10
SEED=42

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
    *) echo "ERROR: unknown REGIME '$REGIME'" >&2; exit 1 ;;
esac

WORK=$CTRL/sims/cov${COV}_n${N_INDIV}_g${N_GEN}_s${SEED}_${SUBDIR_TAG}
READS_DIR=$WORK/reads
for f in $READS_DIR/r1.fq $READS_DIR/r2.fq $WORK/recomb_truth.tsv.gz; do
    [ -s "$f" ] || { echo "ERROR: missing $f -- run 06_run_sim_p80.sh first" >&2; exit 1; }
done

CN_KMER_PREFIX=$CTRL/data/kmer_pa_p80_filt2/kmer_pa
CN_VAR=$CTRL/data/var_pa_p80.var_pa.npz
CN_VAR_META=$CTRL/data/var_pa_p80.meta.npz

OUT_DIR=$CTRL/results/${ODIR}/${REGIME}
mkdir -p $OUT_DIR
SAMPLE=p80_${WTAG}_${REGIME}_cov${COV}_s${SEED}
OUT_TSV=$OUT_DIR/${SAMPLE}.tsv

echo "[$(date)] kMate filt2 GLOBAL weight=$WEIGHT --regime $REGIME"
echo "  reads:   $READS_DIR/r1.fq + r2.fq"
echo "  kmer_pa: $CN_KMER_PREFIX (filt2)"
echo "  weight:  $WEIGHT"
echo "  out:     $OUT_TSV"

$PYTHON -u $DRIVER \
    --kmer-pa-prefix $CN_KMER_PREFIX \
    --var-pa $CN_VAR \
    --var-meta $CN_VAR_META \
    --reads $READS_DIR/r1.fq $READS_DIR/r2.fq \
    --sample $SAMPLE \
    --out $OUT_TSV \
    --threads 8 \
    --chroms Chr1 \
    --block-mode global \
    --kmer-weight $WEIGHT

echo
echo "[$(date)] DONE -- $OUT_TSV"
ls -lh $OUT_TSV
