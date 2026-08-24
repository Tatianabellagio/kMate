#!/bin/bash
#SBATCH --job-name=gen_mat
#SBATCH --account=co_moilab
#SBATCH --partition=savio3_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --output=logs/gen_mat_%A_%a.out
#SBATCH --error=logs/gen_mat_%A_%a.err

# Rebuild analysis/grenenet_gea/gen_matrices/ from the on-cluster af_store.
# 9 array tasks = 3 generations x 3 kinds (snp, nonsnp, smallindel); each task
# writes ONE gen{g}_{kind}_af.npy. Idempotent: a task whose output already
# exists with the right shape is skipped, so preempted tasks resume cleanly.
#
# Partition note: savio4_htc (moilab condo) was DOWN/maintenance at build time,
# so this targets savio3_htc at savio_lowprio (preemptible; --requeue resumes).
#
# Usage:
#   sbatch --array=1-9%9 analysis/grenenet_gea/common/build_gen_matrices_array.sh
#
# Task -> (gen, kind) mapping (1-based array id):
#   id:  1   2     3           4   5     6           7   8     9
#   gen: 1   1     1           2   2     2           3   3     3
#   kind snp nonsnp smallindel snp nonsnp smallindel snp nonsnp smallindel

set -uo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs analysis/grenenet_gea/gen_matrices
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python

GENS=(1 2 3)
KINDS=(snp nonsnp smallindel)
i=$(( SLURM_ARRAY_TASK_ID - 1 ))
GEN=${GENS[$(( i / 3 ))]}
KIND=${KINDS[$(( i % 3 ))]}

echo "[$(date)] task=${SLURM_ARRAY_TASK_ID} gen=${GEN} kind=${KIND}"
$PY -u analysis/grenenet_gea/common/build_gen_matrices.py --gens "$GEN" --kinds "$KIND"
echo "[$(date)] task=${SLURM_ARRAY_TASK_ID} done"
