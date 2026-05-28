#!/bin/bash
# =============================================================================
# launch_loo_downloads.sh
# Submit a SLURM array of LOO downloads. By default downloads ALL 78 manifest
# rows; pass an array spec as $1 to restrict (e.g. "1-5" for the smallest
# subset, "1,3,7" for specific indices).
# =============================================================================
set -euo pipefail
BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping

N=$(awk 'NR>1' $BASE/data/loo_ena_manifest.tsv | wc -l)
ARRAY="${1:-1-${N}%4}"   # default: all, max 4 concurrent (ENA-friendly)

JOB=$(sbatch --parsable --array=$ARRAY $BASE/scripts/download_loo_one.sh)
echo "Submitted LOO download array $JOB ($ARRAY)"
echo ""
echo "After downloads land:"
echo "  bash $BASE/scripts/launch_loo_preprocess.sh '$ARRAY'"
