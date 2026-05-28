#!/bin/bash
#SBATCH --job-name=loo_prep
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=5
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --output=logs/loo_prep_%A_%a.out
#SBATCH --error=logs/loo_prep_%A_%a.err

# =============================================================================
# preprocess_loo_one.sh
# Trimmomatic PE + Clumpify dedup for one LOO sample. Threshold-matched to
# xwu's 1001G individual-accession pipeline at
# TODO: preprocessing history file, not migrating — see PIPELINE_FASTQ_PREPROCESSING.md
# /carnegie/nobackup/scratch/xwu/GrENE_net/vcf/sra/commands.sh:
#   - TruSeq3-PE-2.fa adapter
#   - ILLUMINACLIP 2:30:10:2:True (minAdapterLength=2, keepBothReads=True)
#   - LEADING:5 TRAILING:5 MINLEN:36, no SLIDINGWINDOW
# (xwu's 1001G command is intentionally lighter than the GrENE Pool-seq trim,
# which uses SLIDINGWINDOW:4:20 — older 1001G data was already lower quality
# at the ends and over-trimming there loses too many reads.)
#
# Dedup: Clumpify (BBTools) — fastq-stage equivalent of xwu's Picard
# MarkDuplicates step. We don't have aligned BAMs in this pipeline (PanGenie
# is k-mer based), so a sequence-identity dedup is the appropriate analog.
#
# Reads from data/loo_ena_manifest.tsv, writes to data/loo_preprocessed/.
# =============================================================================
mkdir -p logs
set -eo pipefail
# Activate conda *before* `set -u`: the `pang` env's activation hook
# (cactus_env_vars.sh) does `export PYTHONPATH=...:$PYTHONPATH`, which trips
# nounset when PYTHONPATH is unset on entry.
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pang
set -u

BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping
MANIFEST=$BASE/data/loo_ena_manifest.tsv
RAW_DIR=$BASE/data/loo_raw_fastqs
PREP_DIR=$BASE/data/loo_preprocessed
mkdir -p $PREP_DIR

ADAPT_PE=/global/home/users/tbellg/miniforge3/envs/pang/share/trimmomatic-0.40-0/adapters/TruSeq3-PE-2.fa
ADAPT_SE=/global/home/users/tbellg/miniforge3/envs/pang/share/trimmomatic-0.40-0/adapters/TruSeq3-SE.fa

IDX=${SLURM_ARRAY_TASK_ID:-1}
LINE=$(sed -n "$((IDX+1))p" $MANIFEST)
[ -n "$LINE" ] || { echo "ERROR: no row at $IDX" >&2; exit 1; }
ECOTYPE=$(echo "$LINE" | cut -f1)
RUN=$(echo "$LINE" | cut -f2)

ECO_RAW=$RAW_DIR/$ECOTYPE
RAW_R1=""; RAW_R2=""
if [ -f "$ECO_RAW/${RUN}_1.fastq.gz" ] && [ -f "$ECO_RAW/${RUN}_2.fastq.gz" ]; then
    RAW_R1=$ECO_RAW/${RUN}_1.fastq.gz
    RAW_R2=$ECO_RAW/${RUN}_2.fastq.gz
elif [ -f "$ECO_RAW/${RUN}.fastq.gz" ]; then
    RAW_R1=$ECO_RAW/${RUN}.fastq.gz
else
    echo "ERROR: no raw fastq for $ECOTYPE in $ECO_RAW" >&2; exit 1
fi
echo "[$(date)] LOO $ECOTYPE: R1=$(basename $RAW_R1) R2=$(basename ${RAW_R2:-NONE})"

# OUT_R1/OUT_R2 differ between PE (paired _1P/_2P) and SE (single _dedup). Set
# them after layout detection so the idempotency check matches the actual
# output filenames written by Clumpify below.
if [ -n "$RAW_R2" ]; then
    OUT_R1=$PREP_DIR/${ECOTYPE}_1P_dedup.fq.gz
    OUT_R2=$PREP_DIR/${ECOTYPE}_2P_dedup.fq.gz
    [ -f "$OUT_R1" ] && [ -f "$OUT_R2" ] && { echo "[$(date)] $ECOTYPE already preprocessed (PE), skipping"; exit 0; }
else
    OUT_R1=$PREP_DIR/${ECOTYPE}_dedup.fq.gz
    OUT_R2=""
    [ -f "$OUT_R1" ] && { echo "[$(date)] $ECOTYPE already preprocessed (SE), skipping"; exit 0; }
fi

TRIMMED_DIR=$PREP_DIR/${ECOTYPE}_trim
mkdir -p $TRIMMED_DIR
TRIM_R1=$TRIMMED_DIR/${ECOTYPE}_1P.fq.gz
TRIM_R2=$TRIMMED_DIR/${ECOTYPE}_2P.fq.gz
TRIM_U1=$TRIMMED_DIR/${ECOTYPE}_1U.fq.gz
TRIM_U2=$TRIMMED_DIR/${ECOTYPE}_2U.fq.gz

# Step A: Trimmomatic — exact thresholds from xwu's 1001G command
if [ -n "$RAW_R2" ]; then
    trimmomatic PE -phred33 -threads 5 \
      "$RAW_R1" "$RAW_R2" \
      "$TRIM_R1" "$TRIM_U1" "$TRIM_R2" "$TRIM_U2" \
      ILLUMINACLIP:"$ADAPT_PE":2:30:10:2:True \
      LEADING:5 TRAILING:5 MINLEN:36
else
    TRIM_SE=$TRIMMED_DIR/${ECOTYPE}.trim.fq.gz
    trimmomatic SE -phred33 -threads 5 "$RAW_R1" "$TRIM_SE" \
      ILLUMINACLIP:"$ADAPT_SE":2:30:10 \
      LEADING:5 TRAILING:5 MINLEN:36
fi

# Step B: Clumpify dedup (canonical params)
if [ -n "$RAW_R2" ]; then
    clumpify.sh in1="$TRIM_R1" in2="$TRIM_R2" out1="$OUT_R1" out2="$OUT_R2" \
      dedupe=t dupesubs=0 optical=f -Xmx30g
else
    clumpify.sh in="$TRIM_SE" out="$OUT_R1" \
      dedupe=t dupesubs=0 optical=f -Xmx30g
fi

rm -rf $TRIMMED_DIR
echo "[$(date)] LOO $ECOTYPE preprocess DONE"
ls -lh $OUT_R1 ${OUT_R2:-}
