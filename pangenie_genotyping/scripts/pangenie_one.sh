#!/bin/bash
#SBATCH --job-name=pangenie
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --output=logs/pg_%A_%a.out
#SBATCH --error=logs/pg_%A_%a.err

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
mkdir -p logs
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
eval "$(conda shell.bash hook)"
conda activate pangenie

BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping
MANIFEST=$BASE/data/ena_manifest.tsv
PREP_DIR=$BASE/data/preprocessed
GT_DIR=$BASE/data/genotyped
TMP_DIR=$BASE/data/tmp_genotype
mkdir -p $GT_DIR $TMP_DIR

# Pre-built PanGenie graph index (output of build_pangenie_index.sh).
INDEX_PREFIX=$BASE/data/pang_135_pangenie_index

IDX=${SLURM_ARRAY_TASK_ID:-1}
LINE=$(sed -n "$((IDX+1))p" $MANIFEST)
ECOTYPE=$(echo "$LINE" | cut -f2)

PREP_R1=$PREP_DIR/${ECOTYPE}_1.dedup.fq.gz
PREP_R2=$PREP_DIR/${ECOTYPE}_2.dedup.fq.gz
PREP_SE=$PREP_DIR/${ECOTYPE}.dedup.fq.gz
OUT_PREFIX=$GT_DIR/${ECOTYPE}
OUT_VCF=${OUT_PREFIX}_genotyping.vcf
TMP_FQ=$TMP_DIR/${ECOTYPE}.fq

if [ -f "${OUT_VCF}.gz" ]; then
    echo "[$(date)] $ECOTYPE: already genotyped"; exit 0
fi

# PanGenie requires UNCOMPRESSED reads in a single file. Decompress + concat
# PE pairs into one .fq (PanGenie just k-mer counts, doesn't care about pairing).
trap "rm -f $TMP_FQ" EXIT
echo "[$(date)] $ECOTYPE: decompress reads -> $TMP_FQ"
if [ -f "$PREP_R1" ] && [ -f "$PREP_R2" ]; then
    zcat "$PREP_R1" "$PREP_R2" > "$TMP_FQ"
elif [ -f "$PREP_SE" ]; then
    zcat "$PREP_SE" > "$TMP_FQ"
else
    echo "ERROR: no preprocessed reads for $ECOTYPE" >&2; exit 1
fi
ls -lh "$TMP_FQ"

echo "[$(date)] $ECOTYPE: PanGenie genotype against pang_135"
PanGenie -f $INDEX_PREFIX -i $TMP_FQ -o $OUT_PREFIX -s $ECOTYPE -t 8 -j 8

# PanGenie writes <prefix>_genotyping.vcf — compress + index it. The pangenie
# conda env doesn't ship bgzip/tabix, so reach into the pang env directly.
BGZIP=/global/home/users/tbellg/miniforge3/envs/pang/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/pang/bin/tabix
$BGZIP -f $OUT_VCF
$TABIX -p vcf ${OUT_VCF}.gz

echo "[$(date)] $ECOTYPE: DONE"
ls -lh ${OUT_VCF}.gz
