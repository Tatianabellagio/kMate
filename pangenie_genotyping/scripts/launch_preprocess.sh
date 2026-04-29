#!/bin/bash
# =============================================================================
# launch_preprocess.sh
# Submit a SLURM array of 151 preprocessing jobs (Trimmomatic + Clumpify).
# Auto-detects: each ecotype's source via the manifest.
# Idempotent: skips already-done samples.
# =============================================================================
set -euo pipefail

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping
MANIFEST=$BASE/data/ena_manifest.tsv

N=$(awk -F'\t' 'NR>1' "$MANIFEST" | wc -l)
echo "Submitting array of $N preprocess jobs, 8 concurrent"

JOB=$(sbatch --parsable \
    --array=1-${N}%8 \
    "$BASE/scripts/preprocess_one.sh")
echo "Submitted job array $JOB (1..$N, max 8 concurrent)"
echo ""
echo "Watch progress:"
echo "  squeue -u tbellagio | grep preprocess"
echo "  tail -f $BASE/logs/prep_${JOB}_*.out"
echo ""
echo "After all done, launch PanGenie (Stage 3) — but only after pang_69 + indexing finish:"
echo "  bash $BASE/scripts/launch_pangenie.sh"
