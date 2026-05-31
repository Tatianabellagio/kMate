#!/bin/bash
#SBATCH --job-name=filt2inv_v3
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH --output=logs/filt2inv_v3_%j.out
#SBATCH --error=logs/filt2inv_v3_%j.err
#
# Production kmer_pa filter: filt2 + invariant removal, in one pass FROM RAW.
# keep 2 <= ac <= F-1  (drops ac==0 dead, ac==1 private, ac==F invariant).
# Supersedes build_kmer_pa_filt2_v3qc_v3_chr1.sh (which left ac==F in).
# Output: data/kmer_pa_231_v3qc_v3_filt2inv/cn_<CHR>.{kmer_pa,meta}.npz
#
mkdir -p logs
set -euo pipefail
BASE=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python   # old 'hapfm' env is gone
SRC=$BASE/data/kmer_pa_231_v3qc_v3            # RAW N-on matrix (production base)
OUT=$BASE/data/kmer_pa_231_v3qc_v3_filt2inv
mkdir -p "$OUT"

for CHR in "${@:-Chr1}"; do
    IN_PREFIX=$SRC/cn_${CHR}
    OUT_PREFIX=$OUT/cn_${CHR}
    [ -s "${IN_PREFIX}.kmer_pa.npz" ] || { echo "ERROR: missing ${IN_PREFIX}.kmer_pa.npz"; exit 1; }
    [ ! -s "${OUT_PREFIX}.kmer_pa.npz" ] || { echo "$CHR exists, skipping"; continue; }
    echo "[$(date)] $CHR filt2inv (2 <= ac <= F-1)"
    $PY -u $BASE/src/filter_kmer_pa_production.py \
        --in-prefix "$IN_PREFIX" --out-prefix "$OUT_PREFIX" \
        --min-ac 2 --invariant-margin 1
done
echo "[$(date)] DONE"
