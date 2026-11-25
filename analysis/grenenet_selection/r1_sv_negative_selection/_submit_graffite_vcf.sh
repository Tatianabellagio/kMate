#!/bin/bash
#SBATCH --job-name=gt_vcf
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/gt_vcf_%A.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
/global/home/users/tbellg/miniforge3/envs/basic/bin/python _make_graffite_vcf.py
