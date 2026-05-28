#!/bin/bash
#SBATCH --job-name=chr1_sim_skewed
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/chr1_sim_skewed_%j.out
#SBATCH --error=logs/chr1_sim_skewed_%j.err

set -uo pipefail
cd /global/scratch/users/tbellg/hapfire_sv/poolfreq
mkdir -p tests/logs data/sim_chr1_skewed

# Pick 5 specific founders (use the first 5 from the seqfile for reproducibility)
ALL_F=$(awk 'NR>1 {print $1}' /global/scratch/users/tbellg/pang/pang_1001gplus/pang/seqfile.txt)
FIVE=$(echo "$ALL_F" | head -5 | paste -sd,)
WEIGHTS="0.40,0.25,0.15,0.10,0.10"

echo "[$(date)] Simulating Chr1 SKEWED pool: 5 founders at known weights"
echo "  founders: $FIVE"
echo "  weights:  $WEIGHTS"

/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python tests/simulate_pool.py \
    --founders "$FIVE" \
    --weights "$WEIGHTS" \
    --coverage 30 \
    --read-len 150 \
    --chrom Chr1 \
    --out data/sim_chr1_skewed/skewed5

echo ""
echo "[$(date)] Running pipeline on skewed pool"

# Modify SIM_PREFIX in the e2e test
SIM_PREFIX=data/sim_chr1_skewed/skewed5 \
    /global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -c "
import os, sys
os.environ['SIM_PREFIX_OVERRIDE'] = '$(pwd)/data/sim_chr1_skewed/skewed5'
sys.path.insert(0, 'tests')
exec(open('tests/test_simulation_e2e.py').read().replace(
    'SIM_PREFIX = os.path.join(DATA, \"sim_chr1\", \"uniform82\")',
    f'SIM_PREFIX = os.environ[\"SIM_PREFIX_OVERRIDE\"]'
))
"
