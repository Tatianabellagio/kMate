#!/bin/bash
#SBATCH --job-name=cn_var_v3qc_v3
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --output=logs/cn_var_v3qc_v3_%j.out
#SBATCH --error=logs/cn_var_v3qc_v3_%j.err
mkdir -p logs
set -euo pipefail
BASE=/global/scratch/users/tbellg/hapfire_sv
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
VCF=$BASE/pangenie_genotyping/data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz
OUT_PREFIX=$BASE/poolfreq/data/cn_var_231_v3qc_v3
[ -s "$VCF" ] || { echo "ERROR: missing $VCF"; exit 1; }
[ ! -s "${OUT_PREFIX}.cn_var.npz" ] || exit 0
$PY -u $BASE/poolfreq/src/build_cn_var.py --vcf "$VCF" --out "$OUT_PREFIX"
ls -lh ${OUT_PREFIX}.*
