#!/bin/bash
#SBATCH --job-name=chr1_cnvar
#SBATCH --partition=bse
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/arch3/chr1/A5_cnvar_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/arch3/chr1/A5_cnvar_%j.err
set -euo pipefail

# Phase 2 Job A5: build cn_var matrices from the NEW Chr1 merged panel.
# Outputs: cn_var_231_arch3_chr1.{cn_var,cn_var_called,meta}.npz
# These feed cactus_em via per_sample_per_chrom.py.

cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/arch3/chr1
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
BUILD_CN=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/src/build_cn_var.py
VCF=merged_231_chr1_final.vcf.gz
OUT_PREFIX=cn_var_231_arch3_chr1

[ -s "$VCF" ] || { echo "ERROR: missing $VCF (A4 not done)"; exit 1; }

echo "[$(date)] === build_cn_var on Chr1 merged panel ==="
$PY -u $BUILD_CN --vcf "$VCF" --out "$OUT_PREFIX"
echo
ls -lh ${OUT_PREFIX}.*
echo "[$(date)] DONE A5"
