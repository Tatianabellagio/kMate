set -u
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
RS=/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript
cd "$(dirname "$0")"
$RS export_expression.R $($PY -c "import pandas as pd; print(' '.join(pd.read_csv('results/functional_track_candidates.csv').gene.unique()))")
$PY run_eqtl_cis.py
$PY eqtl_tag_snps.py
$PY plot_shortlist.py --table results/_render_new3.csv
$PY plot_atac.py AT4G02820 AT4G11800 AT2G14910
$PY plot_expression.py AT4G02820 AT4G11800 AT2G14910
$PY organize_figures.py
echo RENDER_DONE
