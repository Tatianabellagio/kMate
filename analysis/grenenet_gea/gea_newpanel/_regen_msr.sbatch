#!/bin/bash
#SBATCH --job-name=msr_chain
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/results/grenenet_gea/gea_newpanel/_regen_logs/msr_chain_%j.log
set -euo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd /global/scratch/users/tbellg/kmate
PYK=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PYB=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
G=analysis/grenenet_gea/gea_newpanel
echo "== run_msr_kendall (snp,nonsnp, bio1) =="
$PYK -u $G/run_msr_kendall.py --classes snp,nonsnp --gen 9 --climate bio1
echo "== compare_msr_snp_nonsnp =="
$PYK -u $G/compare_msr_snp_nonsnp.py
echo "== plot_msr (basic env) =="
$PYB -u $G/plot_msr.py || echo "WARN plot_msr failed (figure only)"
echo DONE_MSR_CHAIN
