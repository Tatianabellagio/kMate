#!/bin/bash
#SBATCH --job-name=seedmix_fix
#SBATCH --partition=savio4_htc
#SBATCH --account=co_moilab
#SBATCH --qos=savio_lowprio
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/qc/seedmix_validation/fix_seedmix/seedmix_%j.log

source ~/miniforge3/etc/profile.d/conda.sh
conda activate kmate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# node-local scratch for the transient jellyfish DB (avoid shared-FS I/O per grenenet/README)
export JF_DIR=${SLURM_TMPDIR:-/tmp}/seedmix_jf_$$
mkdir -p "$JF_DIR"
cd /global/scratch/users/tbellg/kmate
python analysis/grenenet_selection/qc/seedmix_validation/fix_seedmix/run_seedmix_fix.py
rm -rf "$JF_DIR"
