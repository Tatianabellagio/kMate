#!/bin/bash
#SBATCH --job-name=p231_smgrid
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=5:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/benchmarks/p231/logs/07e_smgrid_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/benchmarks/p231/logs/07e_smgrid_%j.err
# 07e -- SEEDMIX (real pool) arms for the h-imbalance grid, on the SAME benchmark
# matrices as the sim columns (kmer_pa_p231 raw + kmer_pa_p231_filt2), so the grid
# is apples-to-apples. Count-once: build one Jellyfish DB from the sample reads,
# query it for all three arms.
#
#   raw      = unfiltered kmer_pa_p231     + --kmer-weight uniform
#   filt2u   = kmer_pa_p231_filt2          + --kmer-weight uniform
#   filt2mb  = kmer_pa_p231_filt2          + --kmer-weight inv_mb   (production)
#
# Usage: bash 07e_run_seedmix_grid_arms.sh S1   (sample id, S1..S8)
set -euo pipefail
S=${1:?Usage: S1..S8}
ROOT=/global/scratch/users/tbellg/kmate
CTRL=$ROOT/benchmarks/p231
PYTHON=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
DRIVER=$ROOT/src/per_sample_per_chrom.py

# Reads from the arch3 manifest (trim-only, PCR-free SEEDMIX).
R1=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/${S}-1.1_P.fq.gz
R2=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/${S}-1.2_P.fq.gz
for f in "$R1" "$R2"; do [ -s "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }; done

RAW_PREFIX=$CTRL/data/kmer_pa_p231/kmer_pa
F2_PREFIX=$CTRL/data/kmer_pa_p231_filt2/kmer_pa
CN_VAR=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
CN_VAR_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz

OUTBASE=$CTRL/results/seedmix_grid
mkdir -p "$OUTBASE"

# ---- count-once: one Jellyfish DB from the read pool ----
DB=${TMPDIR:-/tmp}/seedmix_${S}_grid.jf
echo "[$(date)] build_kmer_db $S -> $DB"
$PYTHON - "$R1" "$R2" "$DB" <<'PY'
import sys
sys.path.insert(0, "/global/scratch/users/tbellg/kmate/src")
from kmer_count import build_kmer_db
r1, r2, db = sys.argv[1], sys.argv[2], sys.argv[3]
build_kmer_db([r1, r2], db, k=31, threads=4)
print("DB built:", db)
PY

run_arm () {  # arm_tag  kmer_prefix  weight
    local TAG=$1 PREFIX=$2 WEIGHT=$3
    local OUT_DIR=$OUTBASE/$TAG; mkdir -p "$OUT_DIR"
    local SAMPLE=SEEDMIX_${S}_${TAG}
    echo "[$(date)] arm=$TAG  weight=$WEIGHT  kmer_pa=$PREFIX"
    $PYTHON -u $DRIVER \
        --kmer-pa-prefix "$PREFIX" \
        --var-pa $CN_VAR --var-meta $CN_VAR_META \
        --reads "$R1" "$R2" --kmer-db "$DB" \
        --sample "$SAMPLE" --out "$OUT_DIR/${SAMPLE}.tsv" \
        --threads 4 --chroms Chr1 --block-mode global --kmer-weight $WEIGHT
}

run_arm raw     "$RAW_PREFIX" uniform
run_arm filt2u  "$F2_PREFIX"  uniform
run_arm filt2mb "$F2_PREFIX"  inv_mb

rm -f "$DB"
echo "[$(date)] DONE SEEDMIX $S  -> $OUTBASE/{raw,filt2u,filt2mb}/SEEDMIX_${S}_*.h_per_chrom.npz"
