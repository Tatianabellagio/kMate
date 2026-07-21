#!/bin/bash
#SBATCH --job-name=nb_binomkend
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/nb_binomkend_%j.log
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/nb_binomkend_%j.log
# Build (kmate) + execute (basic) the SNP-vs-nonSNP binomial+Kendall raw notebook.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel
PYK=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
NBX=/global/home/users/tbellg/miniforge3/envs/basic/bin/jupyter
$PYK -u _build_binom_kendall_nb.py
$NBX nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=7200 notebooks/binom_kendall_snp_vs_nonsnp.ipynb
echo DONE_NB_BINOM_KENDALL
