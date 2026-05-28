#!/bin/bash
# =============================================================================
# launch_pangenie.sh
# Submit a SLURM array of 151 PanGenie genotype jobs.
# Prerequisites:
#   1. preprocess_one.sh has run on all 151 samples (preprocessed/ dir populated)
#   2. build_pangenie_index.sh has run (pang_69 graph indexed)
# =============================================================================
set -euo pipefail

BASE=/global/scratch/users/tbellg/hapfire_sv/pangenie_genotyping
MANIFEST=$BASE/data/ena_manifest.tsv
INDEX=$BASE/data/pang69_pangenie_index

# Sanity checks
if [ ! -f ${INDEX}.cereal ] && [ ! -d $INDEX ]; then
    echo "ERROR: PanGenie index not found ($INDEX). Run build_pangenie_index.sh first." >&2
    exit 1
fi
if [ ! -d $BASE/data/preprocessed ] || [ "$(ls -A $BASE/data/preprocessed 2>/dev/null | wc -l)" -lt 100 ]; then
    echo "ERROR: preprocessed/ dir is empty or sparse. Run preprocess_one.sh array first." >&2
    exit 1
fi

N=$(awk -F'\t' 'NR>1' "$MANIFEST" | wc -l)
echo "Submitting PanGenie array of $N jobs, 8 concurrent"

JOB=$(sbatch --parsable \
    --array=1-${N}%8 \
    "$BASE/scripts/pangenie_one.sh")
echo "Submitted job array $JOB (1..$N, max 8 concurrent)"
echo ""
echo "After all done:"
echo "  bash $BASE/scripts/merge_vcfs.sh   # → 231-founder catalog"
