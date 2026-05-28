#!/bin/bash
#SBATCH --job-name=pg_index
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/pg_index_%j.out
#SBATCH --error=logs/pg_index_%j.err

# =============================================================================
# build_pangenie_index.sh
# Build the PanGenie graph index from pang_69's vcfbub-filtered VCF + ref.
# Run once after pang_69 (job 56180) finishes.
# =============================================================================
mkdir -p logs
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
eval "$(conda shell.bash hook)"
conda activate pangenie

PANG_DIR=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output
PANG_VCF_GZ=$PANG_DIR/pang_1001gplus_all.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa

INDEX_DIR=/global/scratch/users/tbellg/kmate/pangenie_genotyping/data
OUT_PREFIX=$INDEX_DIR/pang_135_pangenie_index
PANG_VCF=$INDEX_DIR/pang_1001gplus_all.dipl.vcf
mkdir -p $INDEX_DIR

if [ ! -s $PANG_VCF_GZ ]; then
    echo "ERROR: cactus VCF not found ($PANG_VCF_GZ) — wait for cactus_all (56180)" >&2
    exit 1
fi

# PanGenie wants (a) uncompressed and (b) diploid genotypes. Cactus pangenome
# writes one haploid GT per founder assembly (`0`, `1`, `.`), but PanGenie's
# GraphBuilder aborts with 'Found invalid genotype. Genotypes must be diploid'.
# A. thaliana ecotypes are inbreds, so each haploid assembly represents a
# homozygous diploid; duplicate the GT (`X` -> `X|X`, `.` -> `.|.`) inline.
# Idempotent: skips if the diploid VCF already exists.
if [ ! -s $PANG_VCF ]; then
    echo "[$(date)] decompress + diploidize $PANG_VCF_GZ -> $PANG_VCF"
    zcat $PANG_VCF_GZ | awk 'BEGIN{OFS="\t"} /^#/{print; next} {for(i=10;i<=NF;i++) $i=$i"|"$i; print}' > $PANG_VCF
fi

echo "[$(date)] PanGenie-index on pang_135 ($(wc -l < $PANG_VCF) lines)"
PanGenie-index -r $REF -v $PANG_VCF -o $OUT_PREFIX -t 8 -k 31

echo "[$(date)] DONE"
ls -lh ${OUT_PREFIX}*
