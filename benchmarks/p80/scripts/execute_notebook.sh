#!/bin/bash
#SBATCH --job-name=p80_nb_exec
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=1:00:00
#SBATCH --output=logs/nb_exec_%j.out
#SBATCH --error=logs/nb_exec_%j.err

# Execute FINAL_RESULTS_cov10_p80.ipynb on a compute node (head node has 16GB
# ulimit which OOM-kills the kernel on the big long-form DataFrame).

mkdir -p logs
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/benchmarks/p80/results

# Notebook tooling lives in the `basic` env (numpy/pandas/scipy/matplotlib +
# jupyter/nbformat/nbconvert/ipykernel). The former `hapfm` env was removed in the
# cluster migration.
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
JUP=/global/home/users/tbellg/miniforge3/envs/basic/bin/jupyter

echo "[$(date)] rebuilding notebook from _build_notebook.py"
$PY _build_notebook.py

echo
echo "[$(date)] executing notebook in-place (timeout 1500s)"
$JUP nbconvert --to notebook --execute --inplace FINAL_RESULTS_cov10_p80.ipynb \
    --ExecutePreprocessor.timeout=2400 \
    --ExecutePreprocessor.kernel_name=python3 \
    --ExecutePreprocessor.startup_timeout=180

echo
echo "[$(date)] DONE"
ls -lh FINAL_RESULTS_cov10_p80.ipynb
