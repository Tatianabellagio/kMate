set -u
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
RS=/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript
cd "$(dirname "$0")"
$RS export_expression.R $($PY -c "import pandas as pd; print(' '.join(pd.read_csv('results/functional_track_candidates.csv').gene.unique()))")
$PY run_eqtl_cis.py
$PY eqtl_tag_snps.py
$PY plot_shortlist.py --table results/_render_new5.csv
$PY plot_atac.py FH4 FRL2 4CL1 AT1G53050 AT1G58235 ZW18 AT3G15570 AT3G44010 AT3G58130 EHD2 BUP HPAT3 AAO AT5G57270 SAC51
$PY plot_expression.py FH4 FRL2 4CL1 AT1G53050 AT1G58235 ZW18 AT3G15570 AT3G44010 AT3G58130 EHD2 BUP HPAT3 AAO AT5G57270 SAC51
$PY organize_figures.py
echo RENDER_DONE
