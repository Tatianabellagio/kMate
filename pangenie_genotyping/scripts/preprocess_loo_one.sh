#!/bin/bash
#SBATCH --job-name=loo_prep
#SBATCH --partition=bse
#SBATCH --cpus-per-task=5
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/loo_prep_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/loo_prep_%A_%a.err

# =============================================================================
# preprocess_loo_one.sh
# Trimmomatic PE + Clumpify dedup for one LOO sample. Conventions match the
# canonical GrENE-Net pipeline at /home/tbellagio/scratch/pang/grenenet_reads/
#   - single `pang` conda env (provides both trimmomatic and clumpify)
#   - same Trimmomatic params (ILLUMINACLIP TruSeq3-PE-2 / SLIDINGWINDOW etc.)
#   - clumpify -Xmx30g, dedupe=t dupesubs=0 optical=f
# Reads from data/loo_ena_manifest.tsv, writes to data/loo_preprocessed/.
# =============================================================================
set -euo pipefail
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pang

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping
MANIFEST=$BASE/data/loo_ena_manifest.tsv
RAW_DIR=$BASE/data/loo_raw_fastqs
PREP_DIR=$BASE/data/loo_preprocessed
mkdir -p $PREP_DIR

ADAPT_PE=/home/tbellagio/miniforge3/envs/pang/share/trimmomatic-0.40-0/adapters/TruSeq3-PE-2.fa
ADAPT_SE=/home/tbellagio/miniforge3/envs/pang/share/trimmomatic-0.40-0/adapters/TruSeq3-SE.fa

IDX=${SLURM_ARRAY_TASK_ID:-1}
LINE=$(sed -n "$((IDX+1))p" $MANIFEST)
[ -n "$LINE" ] || { echo "ERROR: no row at $IDX" >&2; exit 1; }
ECOTYPE=$(echo "$LINE" | cut -f1)
RUN=$(echo "$LINE" | cut -f2)

OUT_R1=$PREP_DIR/${ECOTYPE}_1P_dedup.fq.gz
OUT_R2=$PREP_DIR/${ECOTYPE}_2P_dedup.fq.gz
if [ -f "$OUT_R1" ] && [ -f "$OUT_R2" ]; then
    echo "[$(date)] $ECOTYPE already preprocessed, skipping"; exit 0
fi

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

TRIMMED_DIR=$PREP_DIR/${ECOTYPE}_trim
mkdir -p $TRIMMED_DIR
TRIM_R1=$TRIMMED_DIR/${ECOTYPE}_1P.fq.gz
TRIM_R2=$TRIMMED_DIR/${ECOTYPE}_2P.fq.gz
TRIM_U1=$TRIMMED_DIR/${ECOTYPE}_1U.fq.gz
TRIM_U2=$TRIMMED_DIR/${ECOTYPE}_2U.fq.gz

# Step A: Trimmomatic
if [ -n "$RAW_R2" ]; then
    trimmomatic PE -phred33 -threads 5 \
      "$RAW_R1" "$RAW_R2" \
      "$TRIM_R1" "$TRIM_U1" "$TRIM_R2" "$TRIM_U2" \
      ILLUMINACLIP:"$ADAPT_PE":2:30:10:8:TRUE \
      SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:36
else
    TRIM_SE=$TRIMMED_DIR/${ECOTYPE}.trim.fq.gz
    trimmomatic SE -phred33 -threads 5 "$RAW_R1" "$TRIM_SE" \
      ILLUMINACLIP:"$ADAPT_SE":2:30:10 \
      SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:36
fi

# Step B: Clumpify dedup (canonical params)
if [ -n "$RAW_R2" ]; then
    clumpify.sh in1="$TRIM_R1" in2="$TRIM_R2" out1="$OUT_R1" out2="$OUT_R2" \
      dedupe=t dupesubs=0 optical=f -Xmx30g
else
    clumpify.sh in="$TRIM_SE" out=$PREP_DIR/${ECOTYPE}_dedup.fq.gz \
      dedupe=t dupesubs=0 optical=f -Xmx30g
fi

rm -rf $TRIMMED_DIR
echo "[$(date)] LOO $ECOTYPE preprocess DONE"
ls -lh $OUT_R1 ${OUT_R2:-}
