#!/bin/bash
# =============================================================================
# concat_lanes.sh
# For ecotypes that have multiple per-lane fastqs (the relicts 100001 + 100002
# from the Berkeley rsync), concatenate the 3 lanes per read direction into
# a single _1.fastq.gz / _2.fastq.gz pair matching the ENA layout.
#
# Concatenated gzip streams are valid gzip (gzip-cat is lossless), so we just
# `cat`. Same trick the GrENE-Net trimmed_dedup pipeline uses for multi-lane
# samples.
# =============================================================================
set -euo pipefail

# Repo dir for this stage; override $PANGENIE_GT for sbatch spool copies.
BASE="${PANGENIE_GT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"/data/raw_fastqs

for ECOTYPE in 100001 100002; do
    DIR=$BASE/$ECOTYPE
    OUT_R1=$DIR/${ECOTYPE}_1.fastq.gz
    OUT_R2=$DIR/${ECOTYPE}_2.fastq.gz
    if [ -f "$OUT_R1" ] && [ -f "$OUT_R2" ]; then
        echo "[$(date)] $ECOTYPE: lane-concat already done"
        continue
    fi
    R1_LANES=$(ls $DIR/*_R1_*.fastq.gz 2>/dev/null | sort)
    R2_LANES=$(ls $DIR/*_R2_*.fastq.gz 2>/dev/null | sort)
    if [ -z "$R1_LANES" ]; then
        echo "[$(date)] $ECOTYPE: no per-lane R1 files found, skipping"; continue
    fi
    echo "[$(date)] $ECOTYPE: concat $(echo $R1_LANES | wc -w) R1 lanes + $(echo $R2_LANES | wc -w) R2 lanes"
    cat $R1_LANES > $OUT_R1
    cat $R2_LANES > $OUT_R2
    # Verify
    R1_SIZE=$(ls -la $OUT_R1 | awk '{print $5}')
    R2_SIZE=$(ls -la $OUT_R2 | awk '{print $5}')
    EXPECTED_R1=$(ls -la $R1_LANES | awk '{s += $5} END {print s}')
    EXPECTED_R2=$(ls -la $R2_LANES | awk '{s += $5} END {print s}')
    if [ "$R1_SIZE" != "$EXPECTED_R1" ] || [ "$R2_SIZE" != "$EXPECTED_R2" ]; then
        echo "  ERROR: concat size mismatch" >&2; exit 1
    fi
    # Confirm gzip integrity of the concatenated file
    gzip -t $OUT_R1 2>&1 && gzip -t $OUT_R2 2>&1
    echo "  ✓ $(basename $OUT_R1) ($((R1_SIZE/1024/1024)) MB) + $(basename $OUT_R2) ($((R2_SIZE/1024/1024)) MB)"
done
echo "[$(date)] DONE"
