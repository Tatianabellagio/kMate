#!/bin/bash
#SBATCH --job-name=bigld_install
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:30:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/logs/bigld_install_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/logs/bigld_install_%j.out
set -x
export PATH=/global/home/users/tbellg/miniforge3/bin:$PATH
export CONDA_PKGS_DIRS=/tmp/bigld_pkgs_${SLURM_JOB_ID}
mkdir -p "$CONDA_PKGS_DIRS"
mamba env remove -n bigld -y 2>/dev/null || true
rm -rf /global/home/users/tbellg/miniforge3/envs/bigld
mamba create -y -n bigld -c conda-forge -c bioconda bioconductor-gpart r-optparse
echo "=== verifying gpart ==="
/global/home/users/tbellg/miniforge3/envs/bigld/bin/Rscript -e 'suppressMessages(library(gpart)); cat("GPART_OK\n")'
