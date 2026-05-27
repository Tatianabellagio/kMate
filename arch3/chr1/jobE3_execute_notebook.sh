#!/bin/bash
#SBATCH --job-name=chr1_nbexec
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=00:45:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/arch3/chr1/E3_nbexec_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/arch3/chr1/E3_nbexec_%j.err
set -euo pipefail

# Execute AF_TRUTH_VS_ESTIMATE_arch3_chr1.ipynb end-to-end to verify no bugs.
# Memory budget: cell 2 builds dense (231, 1.9M) arrays = ~1.7 GB each → 3.4 GB peak;
# atomized cn_var load is sparse so ~80 MB; should fit easily in 96G.

cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv
/home/tbellagio/miniforge3/envs/hapfm/bin/jupyter nbconvert \
    --to notebook \
    --execute AF_TRUTH_VS_ESTIMATE_arch3_chr1.ipynb \
    --output AF_TRUTH_VS_ESTIMATE_arch3_chr1.ipynb \
    --ExecutePreprocessor.timeout=600

echo
echo "[$(date)] DONE E3 (notebook executed end-to-end)"
