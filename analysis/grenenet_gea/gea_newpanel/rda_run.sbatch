#!/bin/bash
#SBATCH --job-name=rda_gea
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --account=co_moilab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/results/grenenet_gea/gea_newpanel/rda/slurm_%j.out

set -euo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
DIR=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel

export RDA_CLASSES=snp,nonsnp,sv
echo "===== pRDA (Mahalanobis) ====="
RDA_KSTRUCT=0,2,3 $PY $DIR/rda_prda.py
echo "===== reduced-PC regression (calibration sweep) ====="
RDA_KSTRUCT=0,1,2,3,4,5 $PY $DIR/rda_pcreg.py
echo "DONE"
