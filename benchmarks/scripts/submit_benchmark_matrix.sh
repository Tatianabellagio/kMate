#!/bin/bash
# Submit one run_benchmark_pool job (kMate global + block) per runnable sim pool,
# both panels. Skips the pool given as $1 (already done). Prints submitted job ids.
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
SKIP=${1:-}
for P in p80 p231; do
  for d in benchmarks/$P/sims/*/; do
    pool=$(basename "$d")
    [ -f "$d/reads/r1.fq" ] || continue
    [ -f "$d/recomb_truth.tsv.gz" ] || [ -f "$d/recomb_truth_raw.tsv.gz" ] || continue
    [ "$pool" = "$SKIP" ] && { echo "skip (done): $pool"; continue; }
    jid=$(sbatch --parsable benchmarks/scripts/run_benchmark_pool.sbatch "$pool" "$P")
    echo "$jid  $P  $pool"
  done
done