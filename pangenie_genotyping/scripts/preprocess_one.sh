#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/prep_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/prep_%A_%a.err

# =============================================================================
# preprocess_one.sh
# Stage 2: trim adapters + dedup for one ecotype's raw fastqs.
# Same processing for both ENA and xwu_BAM-derived inputs so all 151 ecotypes
# enter PanGenie at the same logical stage (trimmed + dedup'd paired reads).
#
# Stages:
#   raw fastq → Trimmomatic (adapter+quality trim) → Clumpify (dedup) → ready
#
# Idempotent: skips if output already exists.
#
# Indexed by SLURM_ARRAY_TASK_ID over manifest rows (any source_type).
# =============================================================================
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping
MANIFEST=$BASE/data/ena_manifest.tsv
RAW_DIR=$BASE/data/raw_fastqs
PREP_DIR=$BASE/data/preprocessed
mkdir -p $PREP_DIR

IDX=${SLURM_ARRAY_TASK_ID:-1}
LINE=$(sed -n "$((IDX+1))p" $MANIFEST)
if [ -z "$LINE" ]; then echo "ERROR: no row at $IDX" >&2; exit 1; fi
SRC=$(echo "$LINE" | cut -f1)         # ENA or xwu_BAM
ECOTYPE=$(echo "$LINE" | cut -f2)
RUN=$(echo "$LINE" | cut -f3)
FASTQ_FTP=$(echo "$LINE" | cut -f6)

OUT_R1=$PREP_DIR/${ECOTYPE}_1.dedup.fq.gz
OUT_R2=$PREP_DIR/${ECOTYPE}_2.dedup.fq.gz
if [ -f "$OUT_R1" ] && [ -f "$OUT_R2" ]; then
    echo "[$(date)] $ECOTYPE already preprocessed, skipping"; exit 0
fi

# Step A: get the raw R1/R2 fastqs into a deterministic location
ECO_RAW=$RAW_DIR/$ECOTYPE
mkdir -p $ECO_RAW
RAW_R1=""; RAW_R2=""

if [ "$SRC" = "ENA" ]; then
    # Download script wrote ${RUN}_1.fastq.gz + ${RUN}_2.fastq.gz (or single ${RUN}.fastq.gz)
    if [ -f "$ECO_RAW/${RUN}_1.fastq.gz" ] && [ -f "$ECO_RAW/${RUN}_2.fastq.gz" ]; then
        RAW_R1=$ECO_RAW/${RUN}_1.fastq.gz
        RAW_R2=$ECO_RAW/${RUN}_2.fastq.gz
    elif [ -f "$ECO_RAW/${RUN}.fastq.gz" ]; then
        # Single interleaved fastq — Trimmomatic SE mode, then we use it as a single-end
        # k-mer source (PanGenie / jellyfish handle it fine).
        RAW_R1=$ECO_RAW/${RUN}.fastq.gz
    else
        echo "ERROR: no raw fastq found for $ECOTYPE in $ECO_RAW" >&2; exit 1
    fi
elif [ "$SRC" = "xwu_BAM" ]; then
    # Special case: 100001 / 100002 — use the lane-concat'd outputs from concat_lanes.sh
    RAW_R1=$ECO_RAW/${ECOTYPE}_1.fastq.gz
    RAW_R2=$ECO_RAW/${ECOTYPE}_2.fastq.gz
    if [ ! -f "$RAW_R1" ] || [ ! -f "$RAW_R2" ]; then
        echo "ERROR: lane-concat outputs missing for $ECOTYPE — run concat_lanes.sh first" >&2; exit 1
    fi
else
    echo "ERROR: unknown source_type '$SRC' for $ECOTYPE" >&2; exit 1
fi

echo "[$(date)] $ECOTYPE: source=$SRC R1=$(basename $RAW_R1) R2=$(basename ${RAW_R2:-NONE})"

# Step B: Trimmomatic — same params as GrENE-Net pipeline
#   ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:TRUE  (paired) or
#   ILLUMINACLIP:TruSeq3-SE.fa:2:30:10            (single)
#   SLIDINGWINDOW:4:20  LEADING:5  TRAILING:5  MINLEN:36
TRIM=/home/tbellagio/miniforge3/envs/sequencing_pipeline/share/trimmomatic-0.39-2/trimmomatic.jar
ADAPT_PE=/home/tbellagio/miniforge3/envs/sequencing_pipeline/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa
ADAPT_SE=/home/tbellagio/miniforge3/envs/sequencing_pipeline/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa
TRIMMED_DIR=$PREP_DIR/${ECOTYPE}_trim
mkdir -p $TRIMMED_DIR
TRIM_R1=$TRIMMED_DIR/${ECOTYPE}_1P.fq.gz
TRIM_R2=$TRIMMED_DIR/${ECOTYPE}_2P.fq.gz
TRIM_U1=$TRIMMED_DIR/${ECOTYPE}_1U.fq.gz
TRIM_U2=$TRIMMED_DIR/${ECOTYPE}_2U.fq.gz

if [ -n "$RAW_R2" ]; then
    java -jar $TRIM PE -threads 4 -phred33 \
        $RAW_R1 $RAW_R2 \
        $TRIM_R1 $TRIM_U1 $TRIM_R2 $TRIM_U2 \
        ILLUMINACLIP:${ADAPT_PE}:2:30:10:8:TRUE \
        SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:36
else
    TRIM_R1=$TRIMMED_DIR/${ECOTYPE}.trim.fq.gz
    java -jar $TRIM SE -threads 4 -phred33 \
        $RAW_R1 $TRIM_R1 \
        ILLUMINACLIP:${ADAPT_SE}:2:30:10 \
        SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:36
fi

# Step C: Clumpify dedup (same as GrENE-Net)
CLUMPIFY=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/clumpify.sh
if [ -n "$RAW_R2" ]; then
    $CLUMPIFY in=$TRIM_R1 in2=$TRIM_R2 \
        out=$OUT_R1 out2=$OUT_R2 \
        dedupe=t dupesubs=0 optical=f
else
    $CLUMPIFY in=$TRIM_R1 out=$PREP_DIR/${ECOTYPE}.dedup.fq.gz \
        dedupe=t dupesubs=0 optical=f
fi

# Cleanup intermediate trimmed files
rm -rf $TRIMMED_DIR

echo "[$(date)] $ECOTYPE preprocess DONE"
ls -lh $OUT_R1 ${OUT_R2:-}
