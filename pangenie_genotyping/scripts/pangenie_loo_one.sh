#!/bin/bash
#SBATCH --job-name=loo_pg
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/loo_pg_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/loo_pg_%A_%a.err

# =============================================================================
# pangenie_loo_one.sh
# Stage 7 (LOO variant): genotype one LOO sample with PanGenie against the
# pang_69 catalog, then compare against cactus truth via loo_concordance.py.
#
# Mirrors pangenie_one.sh but reads from data/loo_ena_manifest.tsv +
# data/loo_preprocessed/ and writes to data/loo_genotyped/. After PanGenie
# finishes, runs loo_concordance.py to produce the per-record GC/nRD report.
#
# Indexed by SLURM_ARRAY_TASK_ID over loo_ena_manifest rows (1..78).
# =============================================================================
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
eval "$(conda shell.bash hook)"
conda activate pangenie

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping
MANIFEST=$BASE/data/loo_ena_manifest.tsv
PREP_DIR=$BASE/data/loo_preprocessed
GT_DIR=$BASE/data/loo_genotyped
CONCORD_DIR=$BASE/data/loo_concordance
mkdir -p $GT_DIR $CONCORD_DIR

# pang_69 graph (built once after cactus_all + build_pangenie_index.sh)
PANG69_DIR=/home/tbellagio/scratch/pang/pang_1001gplus/pang_all/output
PANG69_VCF=$PANG69_DIR/pang_1001gplus_all.vcf.gz
REF=/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa

# Cactus assembly-ID → 1001G ID rename map (col 1 = assembly_id, col 2 = ecotype_id)
SAMPLE_RENAME=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/work/sample_rename.txt

IDX=${SLURM_ARRAY_TASK_ID:-1}
LINE=$(sed -n "$((IDX+1))p" $MANIFEST)
[ -n "$LINE" ] || { echo "ERROR: no row at $IDX" >&2; exit 1; }
ECOTYPE=$(echo "$LINE" | cut -f1)

PREP_R1=$PREP_DIR/${ECOTYPE}_1P_dedup.fq.gz
PREP_R2=$PREP_DIR/${ECOTYPE}_2P_dedup.fq.gz
PREP_SE=$PREP_DIR/${ECOTYPE}_dedup.fq.gz
OUT_VCF=$GT_DIR/${ECOTYPE}.vcf

if [ -f "${OUT_VCF}.gz" ] && [ -f "$CONCORD_DIR/${ECOTYPE}_summary.tsv" ]; then
    echo "[$(date)] $ECOTYPE: already genotyped + concordance done"; exit 0
fi

# Build read input string
if [ -f "$PREP_R1" ] && [ -f "$PREP_R2" ]; then
    READS="$PREP_R1 $PREP_R2"
elif [ -f "$PREP_SE" ]; then
    READS="$PREP_SE"
else
    echo "ERROR: no preprocessed reads for $ECOTYPE" >&2; exit 1
fi

# --- Step A: PanGenie genotype --------------------------------------------
if [ ! -f "${OUT_VCF}.gz" ]; then
    echo "[$(date)] $ECOTYPE: PanGenie genotype against pang_69"
    PanGenie -i $READS -r $REF -v $PANG69_VCF -o $OUT_VCF -s $ECOTYPE -t 8 -j 8
    bgzip -f $OUT_VCF
    tabix -p vcf ${OUT_VCF}.gz
else
    echo "[$(date)] $ECOTYPE: ${OUT_VCF}.gz exists, skipping genotype step"
fi

# --- Step B: concordance vs cactus truth ----------------------------------
# Look up the cactus assembly ID for this ecotype (col1 of rename when col2 == ecotype)
TRUTH_SAMPLE=$(awk -v e="$ECOTYPE" '$2==e {print $1}' $SAMPLE_RENAME | head -1)
if [ -z "$TRUTH_SAMPLE" ]; then
    echo "WARN: $ECOTYPE has no cactus assembly mapping in $SAMPLE_RENAME — skipping concordance" >&2
    exit 0
fi
echo "[$(date)] $ECOTYPE: concordance against cactus truth (assembly=$TRUTH_SAMPLE)"

PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python
$PYTHON $BASE/scripts/loo_concordance.py \
    --pangenie-vcf ${OUT_VCF}.gz \
    --truth-vcf    $PANG69_VCF \
    --sample       $ECOTYPE \
    --truth-sample $TRUTH_SAMPLE \
    --out          $CONCORD_DIR/${ECOTYPE}

echo "[$(date)] $ECOTYPE: DONE"
ls -lh ${OUT_VCF}.gz $CONCORD_DIR/${ECOTYPE}_summary.tsv
