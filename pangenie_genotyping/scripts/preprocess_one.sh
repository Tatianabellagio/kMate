#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=logs/prep_%A_%a.out
#SBATCH --error=logs/prep_%A_%a.err

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
mkdir -p logs
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"

BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping
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

# OUT_R1/OUT_R2 differ between PE (paired _1.dedup/_2.dedup) and SE (single
# .dedup). Set them after layout detection so the idempotency check matches
# the actual output filenames written by Clumpify below.
if [ -n "$RAW_R2" ]; then
    OUT_R1=$PREP_DIR/${ECOTYPE}_1.dedup.fq.gz
    OUT_R2=$PREP_DIR/${ECOTYPE}_2.dedup.fq.gz
    [ -f "$OUT_R1" ] && [ -f "$OUT_R2" ] && { echo "[$(date)] $ECOTYPE already preprocessed (PE), skipping"; exit 0; }
else
    OUT_R1=$PREP_DIR/${ECOTYPE}.dedup.fq.gz
    OUT_R2=""
    [ -f "$OUT_R1" ] && { echo "[$(date)] $ECOTYPE already preprocessed (SE), skipping"; exit 0; }
fi

# Step B: Trimmomatic — same params as GrENE-Net pipeline
#   ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:TRUE  (paired) or
#   ILLUMINACLIP:TruSeq3-SE.fa:2:30:10            (single)
#   SLIDINGWINDOW:4:20  LEADING:5  TRAILING:5  MINLEN:36
TRIM=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/share/trimmomatic-0.39-2/trimmomatic.jar
ADAPT_PE=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa
ADAPT_SE=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa
TRIMMED_DIR=$PREP_DIR/${ECOTYPE}_trim
mkdir -p $TRIMMED_DIR
TRIM_R1=$TRIMMED_DIR/${ECOTYPE}_1P.fq.gz
TRIM_R2=$TRIMMED_DIR/${ECOTYPE}_2P.fq.gz
TRIM_U1=$TRIMMED_DIR/${ECOTYPE}_1U.fq.gz
TRIM_U2=$TRIMMED_DIR/${ECOTYPE}_2U.fq.gz

# Trimmomatic — exact thresholds from xwu's 1001G command at
# TODO: preprocessing history file, not migrating — see docs/PIPELINE_FASTQ_PREPROCESSING.md
# /carnegie/nobackup/scratch/xwu/GrENE_net/vcf/sra/commands.sh
# (lighter trim than the GrENE Pool-seq pipeline: minAdapterLength=2 and no
# SLIDINGWINDOW; 1001G data was already lower quality at the ends and
# over-trimming there loses too many reads.)
if [ -n "$RAW_R2" ]; then
    java -jar $TRIM PE -threads 4 -phred33 \
        $RAW_R1 $RAW_R2 \
        $TRIM_R1 $TRIM_U1 $TRIM_R2 $TRIM_U2 \
        ILLUMINACLIP:${ADAPT_PE}:2:30:10:2:True \
        LEADING:5 TRAILING:5 MINLEN:36
else
    TRIM_R1=$TRIMMED_DIR/${ECOTYPE}.trim.fq.gz
    java -jar $TRIM SE -threads 4 -phred33 \
        $RAW_R1 $TRIM_R1 \
        ILLUMINACLIP:${ADAPT_SE}:2:30:10 \
        LEADING:5 TRAILING:5 MINLEN:36
fi

# Step C: Clumpify dedup (same as GrENE-Net)
# Note: sequencing_pipeline env doesn't have clumpify; use pang env's BBTools install
CLUMPIFY=/global/home/users/tbellg/miniforge3/envs/pang/bin/clumpify.sh
if [ -n "$RAW_R2" ]; then
    $CLUMPIFY in=$TRIM_R1 in2=$TRIM_R2 \
        out=$OUT_R1 out2=$OUT_R2 \
        dedupe=t dupesubs=0 optical=f
else
    $CLUMPIFY in=$TRIM_R1 out=$OUT_R1 \
        dedupe=t dupesubs=0 optical=f
fi

# Cleanup intermediate trimmed files
rm -rf $TRIMMED_DIR

echo "[$(date)] $ECOTYPE preprocess DONE"
ls -lh $OUT_R1 ${OUT_R2:-}
