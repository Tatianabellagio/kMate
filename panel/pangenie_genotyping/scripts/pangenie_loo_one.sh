#!/bin/bash
#SBATCH --job-name=loo_pg
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --output=logs/loo_pg_%A_%a.out
#SBATCH --error=logs/loo_pg_%A_%a.err

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
mkdir -p logs
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
eval "$(conda shell.bash hook)"
conda activate pangenie

# Repo dir for this stage; override $PANGENIE_GT for sbatch spool copies.
BASE="${PANGENIE_GT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
MANIFEST=$BASE/data/loo_ena_manifest.tsv
PREP_DIR=$BASE/data/loo_preprocessed
GT_DIR=$BASE/data/loo_genotyped
CONCORD_DIR=$BASE/data/loo_concordance
TMP_DIR=$BASE/data/tmp_loo_genotype
mkdir -p $GT_DIR $CONCORD_DIR $TMP_DIR

# Pre-built PanGenie graph index (one-time output of build_pangenie_index.sh).
INDEX_PREFIX=$BASE/data/pang_135_pangenie_index
# Diploidized cactus VCF — used as ground truth for concordance. Diploid GTs
# match PanGenie's diploid output so loo_concordance.py's dose comparison is
# meaningful (haploid cactus dose=1 vs diploid PanGenie dose=2 would always
# look discordant for inbred carriers).
TRUTH_VCF=$BASE/data/pang_1001gplus_all.dipl.vcf.gz

# Absolute paths to bgzip/tabix because the pangenie env doesn't ship them.
BGZIP=/global/home/users/tbellg/miniforge3/envs/pang/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/pang/bin/tabix

# Cactus assembly-ID → 1001G ID rename map (col 1 = assembly_id, col 2 = ecotype_id).
# REQUIRED external input (see README Prerequisites); override $SAMPLE_RENAME.
SAMPLE_RENAME="${SAMPLE_RENAME:-$BASE/data/cactus_sample_rename.txt}"

IDX=${SLURM_ARRAY_TASK_ID:-1}
LINE=$(sed -n "$((IDX+1))p" $MANIFEST)
[ -n "$LINE" ] || { echo "ERROR: no row at $IDX" >&2; exit 1; }
ECOTYPE=$(echo "$LINE" | cut -f1)

PREP_R1=$PREP_DIR/${ECOTYPE}_1P_dedup.fq.gz
PREP_R2=$PREP_DIR/${ECOTYPE}_2P_dedup.fq.gz
PREP_SE=$PREP_DIR/${ECOTYPE}_dedup.fq.gz
OUT_PREFIX=$GT_DIR/${ECOTYPE}
OUT_VCF=${OUT_PREFIX}_genotyping.vcf
TMP_FQ=$TMP_DIR/${ECOTYPE}.fq

if [ -f "${OUT_VCF}.gz" ] && [ -f "$CONCORD_DIR/${ECOTYPE}_summary.tsv" ]; then
    echo "[$(date)] $ECOTYPE: already genotyped + concordance done"; exit 0
fi

# --- Step A: PanGenie genotype --------------------------------------------
if [ ! -f "${OUT_VCF}.gz" ]; then
    # PanGenie wants UNCOMPRESSED reads in a single file. Decompress + concat.
    trap "rm -f $TMP_FQ" EXIT
    echo "[$(date)] $ECOTYPE: decompress reads -> $TMP_FQ"
    if [ -f "$PREP_R1" ] && [ -f "$PREP_R2" ]; then
        zcat "$PREP_R1" "$PREP_R2" > "$TMP_FQ"
    elif [ -f "$PREP_SE" ]; then
        zcat "$PREP_SE" > "$TMP_FQ"
    else
        echo "ERROR: no preprocessed reads for $ECOTYPE" >&2; exit 1
    fi

    echo "[$(date)] $ECOTYPE: PanGenie genotype against pang_135"
    PanGenie -f $INDEX_PREFIX -i $TMP_FQ -o $OUT_PREFIX -s $ECOTYPE -t 8 -j 8
    $BGZIP -f $OUT_VCF
    $TABIX -p vcf ${OUT_VCF}.gz
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

PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
$PYTHON $BASE/scripts/loo_concordance.py \
    --pangenie-vcf ${OUT_VCF}.gz \
    --truth-vcf    $TRUTH_VCF \
    --sample       $ECOTYPE \
    --truth-sample $TRUTH_SAMPLE \
    --out          $CONCORD_DIR/${ECOTYPE}

echo "[$(date)] $ECOTYPE: DONE"
ls -lh ${OUT_VCF}.gz $CONCORD_DIR/${ECOTYPE}_summary.tsv
