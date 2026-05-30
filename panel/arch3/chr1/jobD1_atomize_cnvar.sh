#!/bin/bash
#SBATCH --job-name=chr1_atomize
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=logs/D1_atomize_%j.out
#SBATCH --error=logs/D1_atomize_%j.err
mkdir -p logs
set -euo pipefail

# Build atomized var_pa from the merged Arch 3 chr1 biallelic VCF.
# Each output row is a single-base substitution (pos, ref_base, alt_base) with
# carriers UNIONed across all source records that imply it. MNPs and overlapping
# region of INS/DEL contribute; pure INS/DEL beyond the alignment overlap do not.

# Run in this script's directory; override $ARCH3_CHR1_DIR when launching from
# an sbatch spool copy outside the source tree.
cd "${ARCH3_CHR1_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
VCF=merged_231_chr1_final.vcf.gz
OUT_PREFIX=var_pa_231_arch3_chr1_atomized

[ -s "$VCF" ] || { echo "ERROR: missing $VCF"; exit 1; }

echo "[$(date)] === atomize var_pa from $VCF ==="
$PY -u build_var_pa_atomized.py --vcf "$VCF" --out "$OUT_PREFIX"
echo
ls -lh ${OUT_PREFIX}.*
echo "[$(date)] DONE D1"
