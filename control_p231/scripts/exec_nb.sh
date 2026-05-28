#!/bin/bash
#SBATCH --job-name=p231_nb
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/control_p231/logs/nb_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/control_p231/logs/nb_%j.err
# Usage: sbatch exec_nb.sh <notebook_path>
set -euo pipefail
NB=${1:?Usage: sbatch exec_nb.sh <notebook_path>}
PY=/home/tbellagio/miniforge3/envs/hapfm/bin
$PY/jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=2400 \
  "$NB"
echo "[$(date)] DONE -- $NB"; ls -lh "$NB"
