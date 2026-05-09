#!/bin/bash
#SBATCH --job-name=hf_blkidx
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/scripts/logs/hf_blkidx_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/scripts/logs/hf_blkidx_%j.err

# One-time: build per-block ecotype-to-haplotype-index lookup from the
# hapFIRE panel + BigLD partition. Used by derive_hapfire_perblock_h.py.
set -uo pipefail
mkdir -p /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/scripts/logs
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq

/usr/bin/time -v /home/tbellagio/miniforge3/envs/hapfm/bin/python -u \
    scripts/build_hapfire_block_index.py 2>&1
echo "[$(date)] DONE"
ls -lh data/hapfire_block_index.npz
