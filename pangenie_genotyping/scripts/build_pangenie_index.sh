#!/bin/bash
#SBATCH --job-name=pg_index
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/pg_index_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/pg_index_%j.err

# =============================================================================
# build_pangenie_index.sh
# Build the PanGenie graph index from pang_69's vcfbub-filtered VCF + ref.
# Run once after pang_69 (job 56180) finishes.
# =============================================================================
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
eval "$(conda shell.bash hook)"
conda activate pangenie

PANG_DIR=/home/tbellagio/scratch/pang/pang_1001gplus/pang_all/output
PANG_VCF_GZ=$PANG_DIR/pang_1001gplus_all.vcf.gz
REF=/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa

INDEX_DIR=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/data
OUT_PREFIX=$INDEX_DIR/pang_135_pangenie_index
PANG_VCF=$INDEX_DIR/pang_1001gplus_all.vcf
mkdir -p $INDEX_DIR

if [ ! -s $PANG_VCF_GZ ]; then
    echo "ERROR: cactus VCF not found ($PANG_VCF_GZ) — wait for cactus_all (56180)" >&2
    exit 1
fi

# PanGenie-index requires an uncompressed VCF (refuses .vcf.gz). Decompress
# once, alongside the index outputs, and reuse on rerun.
if [ ! -s $PANG_VCF ]; then
    echo "[$(date)] decompressing $PANG_VCF_GZ -> $PANG_VCF"
    zcat $PANG_VCF_GZ > $PANG_VCF
fi

echo "[$(date)] PanGenie-index on pang_135 ($(wc -l < $PANG_VCF) lines)"
PanGenie-index -r $REF -v $PANG_VCF -o $OUT_PREFIX -t 8 -k 31

echo "[$(date)] DONE"
ls -lh ${OUT_PREFIX}*
