#!/bin/bash
#SBATCH --job-name=shape_iv
#SBATCH --partition=savio4_htc
#SBATCH --account=co_moilab
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-5
#SBATCH --output=shape_iv_%A_%a.out

# Replicate-aware (two-stage IV / random-effects) response-function GEA.
# idx -> (class, weight): 0 snp/iv 1 snp/ols 2 indel/iv 3 indel/ols 4 sv/iv 5 sv/ols
cd /global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel
source ~/.bashrc
mamba activate kmate
set -eo pipefail
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8

CLS=(snp snp indel indel sv sv)
WGT=(iv ols iv ols iv ols)
i=$SLURM_ARRAY_TASK_ID
echo "task $i -> ${CLS[$i]} / ${WGT[$i]}  on $(hostname)"
python3 shape_gea_iv.py --cls "${CLS[$i]}" --weight "${WGT[$i]}" --n-perm 1000 --seed 0
echo "DONE_SHAPEIV_${CLS[$i]}_${WGT[$i]}"
