#!/bin/bash
#SBATCH --job-name=nb_ksweep
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/nb_ksweep_%j.log
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/nb_ksweep_%j.log
# Build (kmate) + execute (basic) the RAW-LFMM-per-K Manhattan/QQ notebooks (K=3,4,5,6).
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel
PYK=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
NBX=/global/home/users/tbellg/miniforge3/envs/basic/bin/jupyter
for K in 3 4 5 6; do
  $PYK -u _build_ksweep_nb.py $K
  $NBX nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=7200 notebooks/ksweep_K${K}_manhattan_qq.ipynb
done
echo DONE_NB_KSWEEP
