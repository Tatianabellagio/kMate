#!/bin/bash
#SBATCH --job-name=p80_c_em_filt2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/07b_em_filt2_%j.out
#SBATCH --error=logs/07b_em_filt2_%j.err

# =============================================================================
# Phase C (filt2 A/B) -- run cactus_em on p80 sim regime with FILT2 cn_full.
# Mirrors 07_run_cactus_em_p80.sh but:
#   --cn-kmer-prefix .../cn_full_p80_filt2/cn  (singletons dropped)
#   outputs to results/cactus_em_<METHOD>_filt2/<REGIME>/p80_filt2_*.tsv
# sim reads + cn_var + cn_var_called are unchanged (k-mer side only).
#
# Usage:
#   sbatch 07b_run_cactus_em_p80_filt2.sh REGIME METHOD
# where REGIME = n50_g1 | n50_g3 | n80_g1
#       METHOD = global | star2
# =============================================================================
mkdir -p logs
set -euo pipefail

REGIME=${1:?Usage: sbatch 07b_run_cactus_em_p80_filt2.sh REGIME METHOD}
METHOD=${2:?Usage: sbatch 07b_run_cactus_em_p80_filt2.sh REGIME METHOD}

if [[ "$METHOD" != "global" && "$METHOD" != "star2" ]]; then
    echo "ERROR: METHOD must be 'global' or 'star2'; got '$METHOD'" >&2
    exit 1
fi

CTRL=/global/scratch/users/tbellg/hapfire_sv/control_p80
PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
DRIVER=/global/scratch/users/tbellg/hapfire_sv/poolfreq/src/per_sample_per_chrom.py

# Sim work dir (built by 06_run_sim_p80.sh). Seed 42 by convention.
COV=10
SEED=42

case "$REGIME" in
    n50_g1)  N_INDIV=50; N_GEN=1 ;;
    n50_g3)  N_INDIV=50; N_GEN=3 ;;
    n80_g1)  N_INDIV=80; N_GEN=1 ;;
    *) echo "ERROR: unknown REGIME '$REGIME'" >&2; exit 1 ;;
esac

WORK=$CTRL/sims/cov${COV}_n${N_INDIV}_g${N_GEN}_s${SEED}_hotspots_p80_chr1
READS_DIR=$WORK/reads
for f in $READS_DIR/r1.fq $READS_DIR/r2.fq $WORK/recomb_truth.tsv.gz; do
    [ -s "$f" ] || { echo "ERROR: missing $f -- run 06_run_sim_p80.sh first" >&2; exit 1; }
done

CN_KMER_PREFIX=$CTRL/data/cn_full_p80_filt2/cn
CN_VAR=$CTRL/data/cn_var_p80.cn_var.npz
CN_VAR_META=$CTRL/data/cn_var_p80.meta.npz

OUT_DIR=$CTRL/results/cactus_em_${METHOD}_filt2/${REGIME}
mkdir -p $OUT_DIR

SAMPLE=p80_filt2_${REGIME}_cov${COV}_s${SEED}
OUT_TSV=$OUT_DIR/${SAMPLE}.tsv

if [ "$METHOD" = "global" ]; then
    MODE_ARGS="--block-mode global"
else
    MODE_ARGS="--block-mode window --window-bp 10000 --global-anchor-weight 0.3 --hmm-smooth-passes 5 --hmm-smooth-alpha 0.5"
fi

echo "[$(date)] cactus_em FILT2 --method $METHOD --regime $REGIME"
echo "  reads:   $READS_DIR/r1.fq + r2.fq"
echo "  cn_full: $CN_KMER_PREFIX (filt2)"
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
