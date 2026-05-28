#!/bin/bash
#SBATCH --job-name=v3_on_p80
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/08_v3_on_p80_%j.out
#SBATCH --error=logs/08_v3_on_p80_%j.err

# =============================================================================
# Run cactus_em with the V3 panel artifacts against the P80 simulation reads.
# Purpose: apples-to-apples vs p80 panel on the SAME input data, so we can say
# something defensible about "v3 panel vs p80 panel."
#
# NOTE: this writes results that must be evaluated against p80's recomb_truth
# (NOT v3's). The (chrom,pos,ref,alt) intersection of v3's cn_var vs p80's
# cn_var defines the shared evaluation surface in the notebook.
#
# Usage:
#   sbatch 08_run_v3panel_on_p80reads.sh REGIME METHOD
# =============================================================================
mkdir -p logs
set -euo pipefail

REGIME=${1:?Usage: REGIME METHOD}
METHOD=${2:?Usage: REGIME METHOD}

CTRL=/global/scratch/users/tbellg/kmate/control_p80
ROOT=/global/scratch/users/tbellg/kmate
PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
DRIVER=$ROOT/src/per_sample_per_chrom.py

case "$REGIME" in
    n50_g1) N=50; G=1 ;;
    n50_g3) N=50; G=3 ;;
    *) echo "ERROR: unknown REGIME" >&2; exit 1 ;;
esac

WORK=$CTRL/sims/cov10_n50_g${G}_s42_hotspots_p80_chr1
READS_DIR=$WORK/reads
[ -s "$READS_DIR/r1.fq" ] || { echo "ERROR: missing reads -- run 06_run_sim_p80.sh first" >&2; exit 1; }

# V3 panel pointers (unchanged -- 231-founder production panel)
CN_KMER_PREFIX=$ROOT/data/cn_full_231_v3/cn
CN_VAR=$ROOT/data/cn_var_231_v3.cn_var.npz
CN_VAR_META=$ROOT/data/cn_var_231_v3.meta.npz

OUT_DIR=$CTRL/results/v3_panel_on_p80_reads/cactus_em_${METHOD}/${REGIME}
mkdir -p $OUT_DIR
SAMPLE=v3panel_p80reads_${REGIME}_cov10_s42
OUT_TSV=$OUT_DIR/${SAMPLE}.tsv

if [ "$METHOD" = "global" ]; then
    MODE_ARGS="--block-mode global"
else
    MODE_ARGS="--block-mode window --window-bp 10000 --global-anchor-weight 0.3 --hmm-smooth-passes 5 --hmm-smooth-alpha 0.5"
fi

echo "[$(date)] running v3 panel on p80 reads -- $REGIME $METHOD"
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

echo "[$(date)] DONE -- $OUT_TSV"
ls -lh $OUT_TSV
