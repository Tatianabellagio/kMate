#!/bin/bash
#SBATCH --job-name=v3qc_pg
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/v3qc_pg_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/v3qc_pg_%j.err

# =============================================================================
# qc_pg_v4_filter.sh
# v3qc Phase C: build pangenie_153_qc.vcf.gz
#   1. GQ>=20 mask on pangenie_151 + loo_5772 + loo_9947
#   2. merge into pangenie_153
#   3. fill-tags (AC, AN, AC_Het, MAC)
#   4. V4 filter: drop records where AC>0 AND AC_Het/AC >= 0.01 (xwu match)
# =============================================================================
set -euo pipefail

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping
MERGED_DIR=$BASE/data/merged
LOO_DIR=$BASE/data/loo_genotyped
TMP=$BASE/data/v3qc_tmp
OUT_DIR=$BASE/data/v3qc
mkdir -p $TMP $OUT_DIR

BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/tabix

PG151=$MERGED_DIR/pangenie_151.vcf.gz
PG5772=$LOO_DIR/5772_genotyping.vcf.gz
PG9947=$LOO_DIR/9947_genotyping.vcf.gz

echo "[$(date)] v3qc Phase C — build pangenie_153_qc.vcf.gz"
for v in $PG151 $PG5772 $PG9947; do
    [ -s "$v" ] || { echo "ERROR: missing $v"; exit 1; }
done

# ---- Step 1: GQ>=20 mask on each input ----
echo "[$(date)] Step 1: GQ>=20 mask"
PG151_GQ=$TMP/pangenie_151.gq20.vcf.gz
PG5772_GQ=$TMP/5772.gq20.vcf.gz
PG9947_GQ=$TMP/9947.gq20.vcf.gz

for in_vcf in "$PG151:$PG151_GQ" "$PG5772:$PG5772_GQ" "$PG9947:$PG9947_GQ"; do
    src=$(echo "$in_vcf" | cut -d: -f1)
    dst=$(echo "$in_vcf" | cut -d: -f2)
    if [ ! -s "$dst" ]; then
        echo "  $(basename $src) → $(basename $dst)"
        # setGT -e 'FMT/GQ>=20' = DO NOT operate where GQ>=20 (i.e. only mask GQ<20)
        $BCF +setGT $src -Oz -o $dst --threads 4 -- -t q -n . -e 'FMT/GQ>=20'
        $TABIX -p vcf $dst
    else
        echo "  $(basename $dst) exists, skipping"
    fi
done

# ---- Step 2: merge into pangenie_153 ----
echo "[$(date)] Step 2: merge into pangenie_153"
PG153_RAW=$TMP/pangenie_153_raw.vcf.gz
if [ ! -s "$PG153_RAW" ]; then
    $BCF merge $PG151_GQ $PG5772_GQ $PG9947_GQ --threads 8 -Oz -o $PG153_RAW
    $TABIX -p vcf $PG153_RAW
fi
N=$($BCF query -l $PG153_RAW | wc -l)
echo "  samples: $N (expected 153)"
[ "$N" = "153" ] || { echo "ERROR: got $N, expected 153"; exit 1; }

# ---- Step 3: fill-tags ----
echo "[$(date)] Step 3: fill-tags (AC, AN, AC_Het, MAC)"
PG153_FILLED=$TMP/pangenie_153_filled.vcf.gz
if [ ! -s "$PG153_FILLED" ]; then
    # NOTE: bcftools +fill-tags does NOT support MAC. Compute MAC inline downstream
    # via 'AC<2 || (AN-AC)<2' when needed. We only need AC + AC_Het for V4.
    $BCF +fill-tags $PG153_RAW --threads 8 -Oz -o $PG153_FILLED -- -t AC,AN,AC_Het,F_MISSING
    $TABIX -p vcf $PG153_FILLED
fi

# ---- Step 4: V4 < 0.01 filter ----
echo "[$(date)] Step 4: V4 (AC_Het/AC) < 0.01 filter"
PG153_QC=$OUT_DIR/pangenie_153_qc.vcf.gz
# Sanity: count records before
N_BEFORE=$($BCF index -n $PG153_FILLED)
echo "  records before V4 filter: $N_BEFORE"

# Expression: exclude where AC>0 AND AC_Het/AC >= 0.01.
# (Records with AC==0 are kept; the MAC filter will catch them later.)
$BCF view -e '(INFO/AC>0) && ((INFO/AC_Het*1.0/INFO/AC) >= 0.01)' $PG153_FILLED \
    --threads 8 -Oz -o $PG153_QC
$TABIX -p vcf $PG153_QC

N_AFTER=$($BCF index -n $PG153_QC)
echo "  records after V4 filter: $N_AFTER"
echo "  dropped: $((N_BEFORE - N_AFTER)) ($(awk -v b=$N_BEFORE -v a=$N_AFTER 'BEGIN{printf "%.2f%%", 100*(b-a)/b}'))"

echo "[$(date)] DONE — pangenie_153_qc.vcf.gz built"
ls -lh $PG153_QC
