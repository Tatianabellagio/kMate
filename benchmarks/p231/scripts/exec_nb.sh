#!/bin/bash
#SBATCH --job-name=p231_nb
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=2:00:00
#SBATCH --output=logs/nb_%j.out
#SBATCH --error=logs/nb_%j.err
# Usage: sbatch exec_nb.sh <notebook_path>
mkdir -p logs
set -euo pipefail
NB=${1:?Usage: sbatch exec_nb.sh <notebook_path>}
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin
$PY/jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=2400 \
  "$NB"
echo "[$(date)] DONE -- $NB"; ls -lh "$NB"
