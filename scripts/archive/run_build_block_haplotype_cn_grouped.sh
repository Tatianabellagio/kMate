#!/bin/bash
#SBATCH --job-name=build_bldhap_grp
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=logs/build_chr1_mxdiv%j.out
#SBATCH --error=logs/build_chr1_mxdiv%j.err

# Build per-block haplotype-level cn matrix for Chr1 with PHG-style allele
# grouping (k-mer Jaccard distance threshold). Pass MXDIV via --export, or
# defaults to 0.05. Output written to chr1_full_mxdiv<MXDIV>.npz.
#
# Usage:
#   sbatch --export=ALL,MXDIV=0.05 run_build_block_haplotype_cn_grouped.sh
#   sbatch --export=ALL,MXDIV=0.05,MAX_BLOCKS=50 ...   # for quick smoke
mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate

MXDIV=${MXDIV:-0.0}
MAF=${MAF:-0.0}
MAX_BLOCKS=${MAX_BLOCKS:-0}
TAG_PARTS=()
if (( $(echo "$MAF > 0" | bc -l) )); then
    TAG_PARTS+=("maf$(echo $MAF | tr -d '.')")
fi
if (( $(echo "$MXDIV > 0" | bc -l) )); then
    TAG_PARTS+=("mxdiv$(echo $MXDIV | tr -d '.')")
fi
if [[ ${#TAG_PARTS[@]} -eq 0 ]]; then
    TAG="strict"
else
    TAG=$(IFS=_; echo "${TAG_PARTS[*]}")
fi

PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
EXTRA=""
if [[ "$MAX_BLOCKS" -gt 0 ]]; then
    EXTRA="--max-blocks $MAX_BLOCKS"
    OUT=data/block_haplotype_cn/chr1_first${MAX_BLOCKS}_${TAG}.npz
else
    OUT=data/block_haplotype_cn/chr1_full_${TAG}.npz
fi

MAF_ARGS=""
if (( $(echo "$MAF > 0" | bc -l) )); then
    MAF_ARGS="--snp-maf-threshold $MAF \
        --cn-var data/cn_var_231_v2.cn_var.npz \
        --cn-var-meta data/cn_var_231_v2.meta.npz"
fi

echo "[$(date)] starting build_block_haplotype_cn for Chr1 (maf=$MAF, mxdiv=$MXDIV, max_blocks=$MAX_BLOCKS)"
$PYTHON -u scripts/build_block_haplotype_cn.py \
    --block-index sims/visor_freqk/chr1_only_panel/hapfire_block_index_chr1.npz \
    --fastas-dir sims/visor_freqk/founder_fastas_231 \
    --out $OUT \
    --chrom-filter Chr1 --threads 8 \
    --jaccard-mxdiv $MXDIV \
    $MAF_ARGS \
    $EXTRA

echo "[$(date)] DONE → $OUT"
ls -lh $OUT
