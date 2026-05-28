#!/bin/bash
#SBATCH --job-name=build_cnvar
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=2:00:00
#SBATCH --output=logs/build_cnvar_%j.out
#SBATCH --error=logs/build_cnvar_%j.err

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/hapfire_sv

PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
MAX_BLOCKS=${MAX_BLOCKS:-0}
VARIANT_ANCHORED=${VARIANT_ANCHORED:-0}
ALLELE_SPECIFIC_CC1=${ALLELE_SPECIFIC_CC1:-0}
HAMMING=${HAMMING:-0}

EXTRA=""
SUFFIX="cnvar"
if [[ "$VARIANT_ANCHORED" == "1" ]]; then
    EXTRA="$EXTRA --variant-anchored"
    SUFFIX="${SUFFIX}_va"
fi
if [[ "$ALLELE_SPECIFIC_CC1" == "1" ]]; then
    EXTRA="$EXTRA --allele-specific-cc1"
    SUFFIX="${SUFFIX}_cc1"
fi
if [[ "$HAMMING" -gt 0 ]]; then
    EXTRA="$EXTRA --hamming-threshold $HAMMING"
    SUFFIX="${SUFFIX}_h${HAMMING}"
fi
if [[ "$MAX_BLOCKS" -gt 0 ]]; then
    EXTRA="$EXTRA --max-blocks $MAX_BLOCKS"
    OUT=poolfreq/data/block_haplotype_cn/chr1_first${MAX_BLOCKS}_${SUFFIX}.npz
else
    OUT=poolfreq/data/block_haplotype_cn/chr1_full_${SUFFIX}.npz
fi

echo "[$(date)] starting build_block_haplotype_cn_cnvar (max_blocks=$MAX_BLOCKS)"
$PYTHON -u poolfreq/scripts/build_block_haplotype_cn_cnvar.py \
    --block-index sims/visor_freqk/chr1_only_panel/hapfire_block_index_chr1.npz \
    --cn-var poolfreq/data/cn_var_231_v2.cn_var.npz \
    --cn-var-meta poolfreq/data/cn_var_231_v2.meta.npz \
    --fastas-dir sims/visor_freqk/founder_fastas_231 \
    --out $OUT \
    --chrom-filter Chr1 --threads 8 \
    $EXTRA

echo "[$(date)] DONE → $OUT"
ls -lh $OUT
