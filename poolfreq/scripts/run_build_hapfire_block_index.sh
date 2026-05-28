#!/bin/bash
#SBATCH --job-name=hf_blkidx
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/hf_blkidx_%j.out
#SBATCH --error=logs/hf_blkidx_%j.err

# One-time: build per-block ecotype-to-haplotype-index lookup from the
# hapFIRE panel + BigLD partition. Used by derive_hapfire_perblock_h.py.
set -uo pipefail
mkdir -p /global/scratch/users/tbellg/hapfire_sv/poolfreq/scripts/logs
cd /global/scratch/users/tbellg/hapfire_sv/poolfreq

/usr/bin/time -v /global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -u \
    scripts/build_hapfire_block_index.py 2>&1
echo "[$(date)] DONE"
ls -lh data/hapfire_block_index.npz
