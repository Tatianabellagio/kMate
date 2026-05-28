#!/bin/bash
# =============================================================================
# launch_loo_preprocess.sh
# Submit a SLURM array of LOO preprocess jobs. Same array spec semantics as
# launch_loo_downloads.sh — pass $1 to restrict.
# =============================================================================
set -euo pipefail
BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping

N=$(awk 'NR>1' $BASE/data/loo_ena_manifest.tsv | wc -l)
ARRAY="${1:-1-${N}%8}"

JOB=$(sbatch --parsable --array=$ARRAY $BASE/scripts/preprocess_loo_one.sh)
echo "Submitted LOO preprocess array $JOB ($ARRAY)"
