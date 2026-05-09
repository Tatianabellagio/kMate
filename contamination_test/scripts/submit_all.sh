#!/bin/bash
# Submit all (sample, chrom) hapFIRE jobs as a SLURM array.
#
# Usage:
#   bash submit_all.sh [smoke|all]
#     smoke -> just SEEDMIX_S1 chr5 (smallest VCF, fastest validation)
#     all   -> 8 samples x 5 chroms = 40 jobs
#
# Each sub-job runs run_hapfire_perchrom.sh with the corresponding (sample, chrom).

set -eo pipefail

MODE="${1:-smoke}"
ROOT=/home/tbellagio/scratch/hapfire_sv/contamination_test
RUNNER="${ROOT}/scripts/run_hapfire_perchrom.sh"

mkdir -p "${ROOT}/logs"

case "${MODE}" in
  smoke)
    PAIRS=("SEEDMIX_S1 5")
    ;;
  one_sample)
    # One full sample (S1 across all 5 chroms)
    PAIRS=("SEEDMIX_S1 1" "SEEDMIX_S1 2" "SEEDMIX_S1 3" "SEEDMIX_S1 4" "SEEDMIX_S1 5")
    ;;
  all)
    PAIRS=()
    for S in 1 2 3 4 5 6 7 8; do
      for C in 1 2 3 4 5; do
        PAIRS+=("SEEDMIX_S${S} ${C}")
      done
    done
    ;;
  *)
    echo "unknown mode: ${MODE}"
    exit 1
    ;;
esac

# Per-chrom memory (right-sized from observed peaks: chr1=503G, chr2=370G;
# rest estimated from record-count linear scaling + 20% margin).
declare -A MEM=( [1]=600G [2]=450G [3]=500G [4]=450G [5]=550G )

echo "[$(date)] submitting ${#PAIRS[@]} hapFIRE jobs (mode=${MODE})"
for P in "${PAIRS[@]}"; do
  read -r S C <<<"${P}"
  JOB="hapfire_${S}_chr${C}"
  M="${MEM[$C]}"
  ID=$(sbatch --parsable --job-name="${JOB}" --mem="${M}" "${RUNNER}" "${S}" "${C}")
  echo "  ${ID}  ${JOB}  mem=${M}"
done

echo
echo "[$(date)] tip: refresh the stats TSV anytime with"
echo "  bash ${ROOT}/scripts/collect_stats.sh"
echo "  cat  ${ROOT}/results/job_stats.tsv"
