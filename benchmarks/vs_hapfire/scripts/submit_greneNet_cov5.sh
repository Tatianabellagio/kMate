#!/bin/bash
# hapFIRE greneNet arm of the 5x depth column: full poolsize sweep at cov=5x
# (N in {2,5,20,50,150,231} x seeds 42-46), matching the kMate arch3 cov5 arm
# (p231/scripts/submit_poolsize_depth_cov5.sh). Same design as
# submit_greneNet_fair_grid.sh; sim -> hapFIRE, dependency-chained, idempotent.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs

COV=5
POOL_SIZES=(2 5 20 50 150 231)
SEEDS=(42 43 44 45 46)
RES=../results/greneNet_fair

submit_one() {
    local n=$1 seed=$2
    local out="$RES/n${n}_cov${COV}_s${seed}/greneNet_hapfire_n${n}_cov${COV}_s${seed}_ecotype_frequency.txt"
    if [ -s "$out" ]; then
        echo "  [skip] n=$n cov=$COV s=$seed already done"
        return
    fi
    local sim_job hf_job
    sim_job=$(sbatch --parsable run_sim_greneNet.sh "$n" "$COV" "$seed")
    hf_job=$(sbatch --parsable --dependency=afterok:${sim_job} run_hapfire_greneNet.sh "$n" "$COV" "$seed")
    echo "  n=$n cov=$COV s=$seed -> sim=$sim_job hapfire=$hf_job"
}

for n in "${POOL_SIZES[@]}"; do
    for seed in "${SEEDS[@]}"; do
        submit_one "$n" "$seed"
    done
done
