#!/bin/bash
#SBATCH --job-name=af_agree
#SBATCH --partition=bse
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/af_agree_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/af_agree_%j.err
set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv
/home/tbellagio/miniforge3/envs/hapfm/bin/jupyter nbconvert \
  --to notebook --execute --inplace \
  AF_PANEL_AGREEMENT_v3qc_v3_filt2_S1_chr1.ipynb \
  --ExecutePreprocessor.timeout=1800
echo "[$(date)] DONE"
ls -lh plots/AF_PANEL_AGREEMENT_v3qc_v3_filt2_S1_chr1.png 2>&1
