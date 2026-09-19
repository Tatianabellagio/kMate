set -u
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
RS=/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript
cd "$(dirname "$0")"
$RS export_expression.R $($PY -c "import pandas as pd; print(' '.join(pd.read_csv('results/functional_track_candidates.csv').gene.unique()))")
$PY run_eqtl_cis.py
$PY eqtl_tag_snps.py
$PY plot_atac.py SSL7 SCPL34 AT5G40855 CPK32 AT1G30820 AT2G18630 AT5G08460 SKIP27 ASD1
$PY plot_expression.py CYP28 SSL7 SCPL34 AT5G40855 CPK32 AT1G30820 AT2G18630 AT5G08460 SKIP27 ASD1
$PY organize_figures.py
echo RENDER_DONE
