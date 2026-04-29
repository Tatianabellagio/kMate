#!/bin/bash
#SBATCH --job-name=pangenie
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/pg_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/pg_%A_%a.err

# =============================================================================
# pangenie_one.sh
# Stage 3: PanGenie genotype one preprocessed sample against pang_69.
# Outputs a per-sample VCF with GT + DS + GQ for every variant in the cactus
# pangenome catalog.
#
# Indexed by SLURM_ARRAY_TASK_ID over manifest rows.
#
# Prerequisites: pang_69 finished + PanGenie-indexed (genotype_index.cereal etc.)
# =============================================================================
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
eval "$(conda shell.bash hook)"
conda activate pangenie

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping
MANIFEST=$BASE/data/ena_manifest.tsv
PREP_DIR=$BASE/data/preprocessed
GT_DIR=$BASE/data/genotyped
mkdir -p $GT_DIR

# pang_69 graph index (built once after pang_69 finishes)
PANG69_DIR=/home/tbellagio/scratch/pang/pang_1001gplus/pang_all/output
PANG69_VCF=$PANG69_DIR/pang_1001gplus_all.vcf.gz                # vcfbub-filtered
REF=/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa
INDEX_PREFIX=$BASE/data/pang69_pangenie_index    # built by an upstream step

IDX=${SLURM_ARRAY_TASK_ID:-1}
LINE=$(sed -n "$((IDX+1))p" $MANIFEST)
ECOTYPE=$(echo "$LINE" | cut -f2)

PREP_R1=$PREP_DIR/${ECOTYPE}_1.dedup.fq.gz
PREP_R2=$PREP_DIR/${ECOTYPE}_2.dedup.fq.gz
PREP_SE=$PREP_DIR/${ECOTYPE}.dedup.fq.gz
OUT_VCF=$GT_DIR/${ECOTYPE}.vcf

if [ -f "${OUT_VCF}.gz" ]; then
    echo "[$(date)] $ECOTYPE: already genotyped"; exit 0
fi

# Build read input string
if [ -f "$PREP_R1" ] && [ -f "$PREP_R2" ]; then
    READS="$PREP_R1 $PREP_R2"
elif [ -f "$PREP_SE" ]; then
    READS="$PREP_SE"
else
    echo "ERROR: no preprocessed reads for $ECOTYPE" >&2; exit 1
fi

echo "[$(date)] $ECOTYPE: PanGenie genotype against pang_69"

# PanGenie genotype: takes pre-built index + reads + sample name → VCF
PanGenie -i $READS -r $REF -v $PANG69_VCF -o $OUT_VCF -s $ECOTYPE -t 8 -j 8

# Compress + index the output VCF
bgzip -f $OUT_VCF
tabix -p vcf ${OUT_VCF}.gz

echo "[$(date)] $ECOTYPE: DONE"
ls -lh ${OUT_VCF}.gz
