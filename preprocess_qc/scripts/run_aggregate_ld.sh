#!/bin/bash
#SBATCH --job-name=ld_agg
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/preprocess_qc/logs/ld_agg_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/preprocess_qc/logs/ld_agg_%j.err

set -eo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
echo "[$(date)] aggregating ~25 GB of LD pair TSVs"
$PY preprocess_qc/scripts/aggregate_ld.py
echo "[$(date)] DONE"
ls -lh preprocess_qc/output/ld/ld_summary.tsv
