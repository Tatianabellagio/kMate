#!/bin/bash
#SBATCH --job-name=chr1_cnvar
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=logs/A5_cnvar_%j.out
#SBATCH --error=logs/A5_cnvar_%j.err
mkdir -p logs
set -euo pipefail

# Phase 2 Job A5: build var_pa matrices from the NEW Chr1 merged panel.
# Outputs: var_pa_231_arch3_chr1.{var_pa,var_called,meta}.npz
# These feed cactus_em via per_sample_per_chrom.py.

# Run in this script's directory; override $ARCH3_CHR1_DIR when launching from
# an sbatch spool copy outside the source tree.
cd "${ARCH3_CHR1_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
BUILD_CN=../../../src/build_var_pa.py
VCF=merged_231_chr1_final.vcf.gz
OUT_PREFIX=var_pa_231_arch3_chr1

[ -s "$VCF" ] || { echo "ERROR: missing $VCF (A4 not done)"; exit 1; }

echo "[$(date)] === build_var_pa on Chr1 merged panel ==="
$PY -u $BUILD_CN --vcf "$VCF" --out "$OUT_PREFIX"
echo
ls -lh ${OUT_PREFIX}.*
echo "[$(date)] DONE A5"
