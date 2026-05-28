#!/bin/bash
# =============================================================================
# launch_downloads.sh
# Submit a SLURM array of 149 download jobs, max 8 concurrent (ENA-friendly).
# Each task downloads one ecotype's fastq(s) and verifies MD5.
# =============================================================================
set -euo pipefail

BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping
MANIFEST=$BASE/data/ena_manifest.tsv

# Count ENA rows
N=$(awk -F'\t' 'NR>1 && $1=="ENA"' "$MANIFEST" | wc -l)
echo "Submitting array of $N download jobs, %8 concurrent"

mkdir -p $BASE/logs $BASE/data/raw_fastqs

JOB=$(sbatch --parsable \
    --array=1-${N}%8 \
    "$BASE/scripts/download_one.sh")
echo "Submitted job array $JOB (1..$N, max 8 concurrent)"
echo ""
echo "Watch progress:"
echo "  squeue -u tbellg | grep ena_dl"
echo "  tail -f $BASE/logs/dl_${JOB}_*.out"
echo ""
echo "After all done:"
echo "  ls $BASE/data/raw_fastqs/ | wc -l   # should be 149"
echo "  du -sh $BASE/data/raw_fastqs/        # ~152 GB expected"
