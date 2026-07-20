#!/bin/bash
# Extend the p231 kMate poolsize x depth grid (07h_run_kmate_poolsize_depth_p231.sh)
# to N=231 (whole-panel pool) -- same design as N in {2,5,20,50,150}: depth in
# {1,10} x seed in {42-46}. Idempotent -- skips (cov,seed) combos whose kMate
# output already exists. Mirrors speed_vs_hapfire/scripts/submit_greneNet_fair_grid.sh.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs

N=231
DEPTHS=(1 10)
SEEDS=(42 43 44 45 46)

submit_one() {
    local cov=$1 seed=$2
    local out="../results/kmate_chrom_poolsize_depth/n${N}_cov${cov}_s${seed}/p231_chrom_psd_n${N}_cov${cov}_s${seed}.tsv"
    if [ -s "$out" ]; then
        echo "  [skip] cov=$cov s=$seed already done"
        return
    fi
    local sim_job est_job
    sim_job=$(COV=$cov sbatch --parsable 06_run_sim_p231_covarg.sh "$N" 0 "$seed")
    est_job=$(sbatch --parsable --dependency=afterok:${sim_job} 07h_run_kmate_poolsize_depth_p231.sh "$N" "$cov" "$seed")
    echo "  cov=$cov s=$seed -> sim=$sim_job est=$est_job"
}

for cov in "${DEPTHS[@]}"; do
    for seed in "${SEEDS[@]}"; do
        submit_one "$cov" "$seed"
    done
done
