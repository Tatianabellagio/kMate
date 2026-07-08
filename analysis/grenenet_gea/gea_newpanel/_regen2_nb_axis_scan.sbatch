#!/bin/bash
#SBATCH --job-name=nb_axis_scan
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/nb_axis_scan_%j.log
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/nb_axis_scan_%j.log
# Build (kmate) + execute (basic) the site-level axis-scan Manhattan/QQ notebooks (gif + raw).
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel
PYK=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
NBX=/global/home/users/tbellg/miniforge3/envs/basic/bin/jupyter
$PYK -u _build_axis_scan_nb.py pval_gif axis_scan_manhattan_qq.ipynb
$PYK -u _build_axis_scan_nb.py pval_raw axis_scan_manhattan_qq_RAW.ipynb
for NB in axis_scan_manhattan_qq axis_scan_manhattan_qq_RAW; do
  $NBX nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=7200 notebooks/$NB.ipynb
done
echo DONE_NB_AXIS_SCAN
