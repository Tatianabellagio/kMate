#!/bin/bash
#SBATCH --job-name=nb_snpnonsnp
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/nb_snpnonsnp_%j.log
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/nb_snpnonsnp_%j.log
# Build (kmate) + execute (basic) the SNP-vs-nonSNP LFMM+WZA clq0.9 notebook.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel
PYK=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
NBX=/global/home/users/tbellg/miniforge3/envs/basic/bin/jupyter
$PYK -u _build_snp_vs_nonsnp_nb.py
$NBX nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=7200 notebooks/snp_vs_nonsnp_lfmm.ipynb
echo DONE_NB_SNPNONSNP
