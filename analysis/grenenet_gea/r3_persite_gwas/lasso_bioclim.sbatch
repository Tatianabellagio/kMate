#!/bin/bash
#SBATCH --job-name=lasso_bioclim
#SBATCH --partition=savio4_htc
#SBATCH --account=co_moilab
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=logs/lasso_bioclim.%j.out
#SBATCH --error=logs/lasso_bioclim.%j.out

# Design A step 2: multilocus LASSO of founder origin climate (bio1-19), SNP vs SNP+nonSNP.
# HOLD until the GRM/CV pass (varexp_bioclim) confirms the aggregate null. Heavy: genome-wide
# markers x 19 traits x repeated CV of ElasticNetCV. Tune --snp-cap / --nrep down if it overruns.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/analysis/grenenet_gea
source ~/miniforge3/etc/profile.d/conda.sh
conda activate kmate
hostname
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python lasso_bioclim.py --prune-r2 0.9 --window-kb 50 --snp-cap 150000 \
    --l1-ratio 1.0 --nfold 5 --nrep 5
