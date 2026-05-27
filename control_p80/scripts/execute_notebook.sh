#!/bin/bash
#SBATCH --job-name=p80_nb_exec
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=1:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/control_p80/logs/nb_exec_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/control_p80/logs/nb_exec_%j.err

# Execute FINAL_RESULTS_cov10_p80.ipynb on a compute node (head node has 16GB
# ulimit which OOM-kills the kernel on the big long-form DataFrame).

set -euo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/control_p80/results

PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
JUP=/home/tbellagio/miniforge3/envs/hapfm/bin/jupyter

echo "[$(date)] rebuilding notebook from _build_notebook.py"
$PY _build_notebook.py

echo
echo "[$(date)] executing notebook in-place (timeout 1500s)"
$JUP nbconvert --to notebook --execute --inplace FINAL_RESULTS_cov10_p80.ipynb \
    --ExecutePreprocessor.timeout=1500

echo
echo "[$(date)] DONE"
ls -lh FINAL_RESULTS_cov10_p80.ipynb
