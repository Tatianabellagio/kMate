#!/bin/bash
#SBATCH --job-name=audit_cnvar
#SBATCH --partition=bse
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/audit_cnvar_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/audit_cnvar_%j.err
set -euo pipefail
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
$PY -u /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/scripts/audit_cn_var_v3qc_v3.py
