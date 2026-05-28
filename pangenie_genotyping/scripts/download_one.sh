#!/bin/bash
#SBATCH --job-name=ena_dl
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=logs/dl_%A_%a.out
#SBATCH --error=logs/dl_%A_%a.err

# =============================================================================
# Download one ecotype's fastq(s) from ENA, indexed by SLURM_ARRAY_TASK_ID.
# Verifies MD5 if available. Idempotent: skips if file exists with correct MD5.
#
# The array index 1..N maps to the N-th ENA-source row in the manifest
# (xwu_BAM rows are skipped — those don't need downloading).
# =============================================================================
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping
MANIFEST=$BASE/data/ena_manifest.tsv
OUT_DIR=$BASE/data/raw_fastqs
mkdir -p $OUT_DIR

# Parse ENA-only rows (skip header + xwu_BAM rows). Use the array index 1-based.
IDX=${SLURM_ARRAY_TASK_ID:-1}
LINE=$(awk -F'\t' 'NR==1 || $1=="ENA"' "$MANIFEST" | sed -n "$((IDX+1))p")
if [ -z "$LINE" ]; then
    echo "ERROR: no row at index $IDX" >&2; exit 1
fi

ECOTYPE=$(echo "$LINE" | cut -f2)
RUN=$(echo "$LINE" | cut -f3)
FASTQ_FTP=$(echo "$LINE" | cut -f6)
FASTQ_MD5=$(echo "$LINE" | cut -f7)

echo "[$(date)] ecotype=$ECOTYPE run=$RUN files=$FASTQ_FTP"

ECO_DIR=$OUT_DIR/$ECOTYPE
mkdir -p $ECO_DIR

# Split URLs and MD5s on ; (semi-colon-separated lists)
IFS=';' read -ra URLS <<< "$FASTQ_FTP"
IFS=';' read -ra MD5S <<< "$FASTQ_MD5"

for i in "${!URLS[@]}"; do
    URL="${URLS[$i]}"
    EXP_MD5="${MD5S[$i]:-}"
    # ENA gives the FTP host without scheme — prepend https:// (faster than ftp://
    # over WAN and supported by ENA at https://ftp.sra.ebi.ac.uk/...)
    if [[ ! "$URL" =~ ^(ftp|http)s?:// ]]; then
        URL="https://$URL"
    fi
    BASENAME=$(basename "$URL")
    OUT="$ECO_DIR/$BASENAME"

    # Skip if already correct
    if [ -f "$OUT" ] && [ -n "$EXP_MD5" ]; then
        ACT_MD5=$(md5sum "$OUT" | cut -d' ' -f1)
        if [ "$ACT_MD5" = "$EXP_MD5" ]; then
            echo "  [skip] $OUT (md5 matches)"
            continue
        else
            echo "  [redo] $OUT (md5 mismatch: got $ACT_MD5, want $EXP_MD5)"
            rm -f "$OUT"
        fi
    elif [ -f "$OUT" ]; then
        echo "  [skip] $OUT (file exists, no md5 to verify)"
        continue
    fi

    # Download with wget — simpler, with --continue for resume support.
    # ENA rate-limits, so a single connection is realistic anyway.
    echo "  [dl] $URL"
    wget --quiet --continue --tries=5 --timeout=120 \
        -O "$OUT" \
        "$URL"

    # Verify MD5
    if [ -n "$EXP_MD5" ]; then
        ACT_MD5=$(md5sum "$OUT" | cut -d' ' -f1)
        if [ "$ACT_MD5" != "$EXP_MD5" ]; then
            echo "  ERROR: md5 mismatch on $OUT" >&2
            exit 1
        fi
        echo "  [ok] md5 verified"
    fi
done

# Report disk usage
du -sh "$ECO_DIR"
echo "[$(date)] DONE ecotype=$ECOTYPE"
