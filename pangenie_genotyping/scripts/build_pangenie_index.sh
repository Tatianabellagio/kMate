#!/bin/bash
#SBATCH --job-name=pg_index
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
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

PANG69_DIR=/home/tbellagio/scratch/pang/pang_1001gplus/pang_all/output
PANG69_VCF=$PANG69_DIR/pang_1001gplus_all.vcf.gz
REF=/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa

OUT_PREFIX=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/data/pang69_pangenie_index
mkdir -p $(dirname $OUT_PREFIX)

if [ ! -s $PANG69_VCF ]; then
    echo "ERROR: pang_69 VCF not found ($PANG69_VCF) — wait for cactus_all (56180)" >&2
    exit 1
fi

echo "[$(date)] PanGenie-index pang_69"
PanGenie-index -r $REF -v $PANG69_VCF -o $OUT_PREFIX -t 8 -k 31

echo "[$(date)] DONE"
ls -lh ${OUT_PREFIX}*
