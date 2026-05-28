#!/bin/bash
#SBATCH --job-name=chr1_sim
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/chr1_sim_%j.out
#SBATCH --error=logs/chr1_sim_%j.err

set -uo pipefail
cd /global/scratch/users/tbellg/hapfire_sv/poolfreq
mkdir -p tests/logs data/sim_chr1

# Build comma-separated lists of all 82 founder Assembly_IDs and uniform weights
FOUNDER_LIST=$(awk 'NR>1 {print $1}' /global/scratch/users/tbellg/pang/pang_1001gplus/pang/seqfile.txt | paste -sd,)
N_FOUNDERS=$(awk 'NR>1' /global/scratch/users/tbellg/pang/pang_1001gplus/pang/seqfile.txt | wc -l)
WEIGHTS=$(yes "1" | head -n $N_FOUNDERS | paste -sd,)

echo "[$(date)] Simulating Chr1 pool: $N_FOUNDERS uniform-weight founders at 30× coverage"

/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python tests/simulate_pool.py \
    --founders "$FOUNDER_LIST" \
    --weights "$WEIGHTS" \
    --coverage 30 \
    --read-len 150 \
    --chrom Chr1 \
    --out data/sim_chr1/uniform82

echo ""
echo "[$(date)] Running pipeline on simulated pool"

# Now run our pipeline against this pool
/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python tests/test_simulation_e2e.py
