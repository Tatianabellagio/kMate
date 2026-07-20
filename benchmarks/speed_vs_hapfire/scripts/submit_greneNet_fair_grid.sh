#!/bin/bash
# Submit the fair hapFIRE (greneNet-panel) sim -> hapFIRE chain across the FULL
# poolsize x depth grid, matching kMate's existing arch3 benchmark
# (benchmarks/p231/results/kmate_chrom_poolsize_depth/): N in {2,5,20,50,150,231} at
# cov in {1,10}, plus N=50 at cov in {30,50}, x seeds {42-46}. Idempotent -- 
# skips (N,cov,seed) combos whose hapFIRE output already exists.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs

POOL_SIZES=(2 5 20 50 150 231)
DEPTHS_ALL=(1 10)
DEPTHS_N50_EXTRA=(30 50)
SEEDS=(42 43 44 45 46)
RES=../results/greneNet_fair

submit_one() {
    local n=$1 cov=$2 seed=$3
    local out="$RES/n${n}_cov${cov}_s${seed}/greneNet_hapfire_n${n}_cov${cov}_s${seed}_ecotype_frequency.txt"
    if [ -s "$out" ]; then
        echo "  [skip] n=$n cov=$cov s=$seed already done"
        return
    fi
    local sim_job hf_job
    sim_job=$(sbatch --parsable run_sim_greneNet.sh "$n" "$cov" "$seed")
    hf_job=$(sbatch --parsable --dependency=afterok:${sim_job} run_hapfire_greneNet.sh "$n" "$cov" "$seed")
    echo "  n=$n cov=$cov s=$seed -> sim=$sim_job hapfire=$hf_job"
}

for n in "${POOL_SIZES[@]}"; do
    for cov in "${DEPTHS_ALL[@]}"; do
        for seed in "${SEEDS[@]}"; do
            submit_one "$n" "$cov" "$seed"
        done
    done
done
for cov in "${DEPTHS_N50_EXTRA[@]}"; do
    for seed in "${SEEDS[@]}"; do
        submit_one 50 "$cov" "$seed"
    done
done
