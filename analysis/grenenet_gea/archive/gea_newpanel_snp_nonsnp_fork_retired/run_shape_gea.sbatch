#!/bin/bash
#SBATCH --job-name=shape_gea
#SBATCH --partition=savio4_htc
#SBATCH --account=co_moilab
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --array=0-5
#SBATCH --output=shape_gea_%A_%a.out

# Response-function (shape) GEA over all classes x weightings.
# Array index -> (class, weight): 0 snp/none 1 snp/flowers 2 indel/none
#                                 3 indel/flowers 4 sv/none 5 sv/flowers
cd /global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel
source ~/.bashrc
mamba activate kmate
set -eo pipefail
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8

CLS=(snp snp indel indel sv sv)
WGT=(none flowers none flowers none flowers)
i=$SLURM_ARRAY_TASK_ID
echo "task $i -> ${CLS[$i]} / ${WGT[$i]}  on $(hostname)"
python3 shape_gea.py --cls "${CLS[$i]}" --weight "${WGT[$i]}" --n-perm 1000 --seed 0
echo "DONE_SHAPE_${CLS[$i]}_${WGT[$i]}"
