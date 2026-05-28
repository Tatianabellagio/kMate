#!/bin/bash
#SBATCH --job-name=cn_var_v3qc_v2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --output=logs/cn_var_v3qc_v2_%j.out
#SBATCH --error=logs/cn_var_v3qc_v2_%j.err
mkdir -p logs
set -euo pipefail
BASE=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
VCF=$BASE/panel/pangenie_genotyping/data/v3qc_v2/founders_231_v3qc_v2.haploid.vcf.gz
OUT_PREFIX=$BASE/data/cn_var_231_v3qc_v2
[ -s "$VCF" ] || { echo "ERROR: missing $VCF"; exit 1; }
[ ! -s "${OUT_PREFIX}.cn_var.npz" ] || exit 0
$PY -u $BASE/src/build_cn_var.py --vcf "$VCF" --out "$OUT_PREFIX"
ls -lh ${OUT_PREFIX}.cn_var.npz ${OUT_PREFIX}.meta.npz
