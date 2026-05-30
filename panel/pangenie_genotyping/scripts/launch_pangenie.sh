#!/bin/bash
# =============================================================================
# launch_pangenie.sh
# Submit a SLURM array of 151 PanGenie genotype jobs.
# Prerequisites:
#   1. preprocess_one.sh has run on all 151 samples (preprocessed/ dir populated)
#   2. build_pangenie_index.sh has run (pang_135 graph indexed)
# =============================================================================
set -euo pipefail

# Repo dir for this stage; override $PANGENIE_GT for sbatch spool copies.
BASE="${PANGENIE_GT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
MANIFEST=$BASE/data/ena_manifest.tsv
# Per-chrom prefix produced by build_pangenie_index.sh; PanGenie appends _Chr<N>_*.
# (Consumed only by the preflight check below — pangenie_one.sh hardcodes its own prefix.)
INDEX=$BASE/data/pang_135_pangenie_index

# Sanity checks
if [ ! -f ${INDEX}_Chr1_Graph.cereal ]; then
    echo "ERROR: PanGenie index not found (${INDEX}_Chr1_Graph.cereal). Run build_pangenie_index.sh first." >&2
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
