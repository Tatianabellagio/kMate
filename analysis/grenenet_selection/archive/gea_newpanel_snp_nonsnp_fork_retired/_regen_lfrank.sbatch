#!/bin/bash
#SBATCH --job-name=lfrank_chain
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/results/_regen_logs/lfrank_chain_%j.log
set -euo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd /global/scratch/users/tbellg/kmate
PYK=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PYB=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
G=analysis/grenenet_gea/gea_newpanel
for CLS in snp nonsnp; do
  echo "== run_lf_rank $CLS (gen9 bio1 K16) =="
  $PYK -u $G/run_lf_rank.py --class "$CLS" --gen 9 --climate bio1 --k 16
done
echo "== lf_rank_compare =="
$PYK -u $G/lf_rank_compare.py
echo "== lf_rank_fig (basic env) =="
$PYB -u $G/lf_rank_fig.py || echo "WARN lf_rank_fig failed (figure only)"
echo DONE_LFRANK_CHAIN
