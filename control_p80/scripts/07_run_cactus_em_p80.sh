#!/bin/bash
#SBATCH --job-name=p80_c_em
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/07_em_%j.out
#SBATCH --error=logs/07_em_%j.err

# =============================================================================
# Phase C -- run cactus_em (one method) on one p80 sim regime.
#
# Usage:
#   sbatch 07_run_cactus_em_p80.sh REGIME METHOD
# where REGIME = n50_g1 | n50_g3 | n80_g1
#       METHOD = global | star2
#
# Outputs go to control_p80/results/<METHOD>/<REGIME>/<sample>.tsv
# =============================================================================
mkdir -p logs
set -euo pipefail

REGIME=${1:?Usage: sbatch 07_run_cactus_em_p80.sh REGIME METHOD}
METHOD=${2:?Usage: sbatch 07_run_cactus_em_p80.sh REGIME METHOD}

if [[ "$METHOD" != "global" && "$METHOD" != "star2" ]]; then
    echo "ERROR: METHOD must be 'global' or 'star2'; got '$METHOD'" >&2
    exit 1
fi

CTRL=/global/scratch/users/tbellg/kmate/control_p80
PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
DRIVER=/global/scratch/users/tbellg/kmate/src/per_sample_per_chrom.py

# Sim work dir (built by 06_run_sim_p80.sh). Seed 42 by convention.
COV=10
SEED=42

case "$REGIME" in
    n50_g0)         N_INDIV=50;  N_GEN=0; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n80_g0)         N_INDIV=80;  N_GEN=0; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n200_g0)        N_INDIV=200; N_GEN=0; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n231_g0)        N_INDIV=231; N_GEN=0; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n50_g1)         N_INDIV=50;  N_GEN=1; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n50_g3)         N_INDIV=50;  N_GEN=3; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n80_g1)         N_INDIV=80;  N_GEN=1; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n200_g1)        N_INDIV=200; N_GEN=1; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n231_g1)        N_INDIV=231; N_GEN=1; SUBDIR_TAG="hotspots_p80_chr1" ;;
    n50_g3_dom500)  N_INDIV=50;  N_GEN=3; SUBDIR_TAG="hotspots_dom500_p80_chr1" ;;
    *) echo "ERROR: unknown REGIME '$REGIME'" >&2; exit 1 ;;
esac

WORK=$CTRL/sims/cov${COV}_n${N_INDIV}_g${N_GEN}_s${SEED}_${SUBDIR_TAG}
READS_DIR=$WORK/reads
for f in $READS_DIR/r1.fq $READS_DIR/r2.fq $WORK/recomb_truth.tsv.gz; do
    [ -s "$f" ] || { echo "ERROR: missing $f -- run 06_run_sim_p80.sh first" >&2; exit 1; }
done

CN_KMER_PREFIX=$CTRL/data/cn_full_p80/cn
CN_VAR=$CTRL/data/cn_var_p80.cn_var.npz
CN_VAR_META=$CTRL/data/cn_var_p80.meta.npz

OUT_DIR=$CTRL/results/cactus_em_${METHOD}/${REGIME}
mkdir -p $OUT_DIR

SAMPLE=p80_${REGIME}_cov${COV}_s${SEED}
OUT_TSV=$OUT_DIR/${SAMPLE}.tsv

if [ "$METHOD" = "global" ]; then
    MODE_ARGS="--block-mode global"
else
    MODE_ARGS="--block-mode window --window-bp 10000 --global-anchor-weight 0.3 --hmm-smooth-passes 5 --hmm-smooth-alpha 0.5"
fi

echo "[$(date)] cactus_em --method $METHOD --regime $REGIME"
echo "  reads:   $READS_DIR/r1.fq + r2.fq"
echo "  cn_full: $CN_KMER_PREFIX"
echo "  cn_var:  $CN_VAR"
echo "  out:     $OUT_TSV"
echo "  mode:    $MODE_ARGS"

$PYTHON -u $DRIVER \
    --cn-kmer-prefix $CN_KMER_PREFIX \
    --cn-var $CN_VAR \
    --cn-var-meta $CN_VAR_META \
    --reads $READS_DIR/r1.fq $READS_DIR/r2.fq \
    --sample $SAMPLE \
    --out $OUT_TSV \
    --threads 8 \
    --chroms Chr1 \
    $MODE_ARGS

echo
echo "[$(date)] DONE -- $OUT_TSV"
ls -lh $OUT_TSV
