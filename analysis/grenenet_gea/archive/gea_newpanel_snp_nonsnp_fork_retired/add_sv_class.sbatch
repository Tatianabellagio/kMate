#!/bin/bash
#SBATCH --job-name=sv_axis_scan
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/sv_axis_scan_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/sv_axis_scan_%j.out
# Add the SV-ONLY class to the site-level LFMM axis scan (bio1..19 + pc1), then
# rebuild the 3-class (snp / nonsnp / sv) significant-blocks table + overlaps.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
NP=analysis/grenenet_gea/gea_newpanel
SITE=analysis/grenenet_gea/gea_newpanel/results/lfmm_site
ES=analysis/grenenet_gea/gea_newpanel/results/env_site
POW=analysis/grenenet_gea/gea_newpanel/results/power

echo "== build SV-only site-level Δp matrix =="
$PY -u $NP/build_power_inputs.py --classes sv --skip-env

echo "== SV site-level LFMM scan bio1..19 + pc1 (K=3, raw+GIF) =="
for AX in $(seq -f 'bio%g' 1 19) pc1; do
  E="$ES/env_site_${AX}.csv"; [ -f "$E" ] || continue
  $R $NP/run_lfmm_both_p.R "$SITE/lfmm_sv_site" "$E" 3 "$POW/exp2_scan_sv_${AX}_bothp.csv"
done

echo "== rebuild 3-class significant-blocks table (snp/nonsnp/sv) + overlap =="
$PY -u $NP/sig_blocks_table.py
echo DONE_SV_CLASS
