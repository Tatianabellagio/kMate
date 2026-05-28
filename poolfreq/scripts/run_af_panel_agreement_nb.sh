#!/bin/bash
#SBATCH --job-name=af_agree
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=logs/af_agree_%j.out
#SBATCH --error=logs/af_agree_%j.err
mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/hapfire_sv
/global/home/users/tbellg/miniforge3/envs/hapfm/bin/jupyter nbconvert \
  --to notebook --execute --inplace \
  AF_PANEL_AGREEMENT_v3qc_v3_filt2_S1_chr1.ipynb \
  --ExecutePreprocessor.timeout=1800
echo "[$(date)] DONE"
ls -lh plots/AF_PANEL_AGREEMENT_v3qc_v3_filt2_S1_chr1.png 2>&1
