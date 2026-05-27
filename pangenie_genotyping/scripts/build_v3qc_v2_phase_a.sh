#!/bin/bash
#SBATCH --job-name=v3qc_v2_A
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/v3qc_v2_A_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/v3qc_v2_A_%j.err

# Phase A: rebuild v3qc with decompose-first ordering.
# Order: GQ mask (already done) → norm-decompose → fill-tags → V4 filter.
# No PG-MAC step here — we'll decide that AFTER seeing the post-V4 MAC distribution.
#
# Output:
#   cactus_78_bi.vcf.gz                — biallelic cactus_78 (78 cactus founders)
#   pangenie_153_qc_v2.vcf.gz          — biallelic PG, GQ+V4 only (no MAC)
set -euo pipefail

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/data
BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/tabix
REF=/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa

CACTUS78=$BASE/v3qc_tmp/cactus_78.vcf.gz
PG_RAW=$BASE/v3qc_tmp/pangenie_153_raw.vcf.gz   # post-GQ-mask + post-merge PG, pre-V4
OUT_DIR=$BASE/v3qc_v2
mkdir -p $OUT_DIR

CACTUS78_BI=$OUT_DIR/cactus_78_bi.vcf.gz
PG_BI=$OUT_DIR/pangenie_153_bi.vcf.gz
PG_FILLED=$OUT_DIR/pangenie_153_filled_bi.vcf.gz
PG_QC=$OUT_DIR/pangenie_153_qc_v2.vcf.gz

[ -s "$CACTUS78" ] || { echo "ERROR: missing $CACTUS78"; exit 1; }
[ -s "$PG_RAW" ] || { echo "ERROR: missing $PG_RAW"; exit 1; }
[ -s "$REF" ] || { echo "ERROR: missing $REF"; exit 1; }

# Step 1: decompose cactus_78 to biallelic
echo "[$(date)] Step 1: norm cactus_78 → biallelic"
if [ ! -s "$CACTUS78_BI" ]; then
    $BCF norm -f $REF -m -any $CACTUS78 --threads 8 -Ou 2> $OUT_DIR/cactus_norm.log | \
        $BCF sort -m 8G -Oz -o $CACTUS78_BI -
    $TABIX -p vcf $CACTUS78_BI
fi
N_C_IN=$($BCF index -n $CACTUS78)
N_C_OUT=$($BCF index -n $CACTUS78_BI)
echo "  cactus_78: $N_C_IN → $N_C_OUT after decompose"

# Step 2: decompose pangenie_153 to biallelic
echo "[$(date)] Step 2: norm pangenie_153 → biallelic"
if [ ! -s "$PG_BI" ]; then
    $BCF norm -f $REF -m -any $PG_RAW --threads 8 -Ou 2> $OUT_DIR/pg_norm.log | \
        $BCF sort -m 8G -Oz -o $PG_BI -
    $TABIX -p vcf $PG_BI
fi
N_P_IN=$($BCF index -n $PG_RAW)
N_P_OUT=$($BCF index -n $PG_BI)
echo "  pangenie_153: $N_P_IN → $N_P_OUT after decompose"

# Step 3: fill-tags on biallelic PG
echo "[$(date)] Step 3: fill-tags AC, AN, AC_Het, F_MISSING on biallelic PG"
if [ ! -s "$PG_FILLED" ]; then
    $BCF +fill-tags $PG_BI --threads 8 -Oz -o $PG_FILLED -- -t AC,AN,AC_Het,F_MISSING
    $TABIX -p vcf $PG_FILLED
fi

# Step 4: V4 filter (per-ALT correct now since biallelic)
echo "[$(date)] Step 4: V4 filter (AC_Het/AC < 0.01)"
N_BEFORE=$($BCF index -n $PG_FILLED)
$BCF view -e '(INFO/AC>0) && ((INFO/AC_Het*1.0/INFO/AC) >= 0.01)' $PG_FILLED --threads 8 -Oz -o $PG_QC
$TABIX -p vcf $PG_QC
N_AFTER=$($BCF index -n $PG_QC)
echo "  PG records before V4:  $N_BEFORE"
echo "  PG records after V4:   $N_AFTER"
echo "  Dropped: $((N_BEFORE - N_AFTER)) ($(awk -v b=$N_BEFORE -v a=$N_AFTER 'BEGIN{printf "%.2f%%", 100*(b-a)/b}'))"

echo ""
echo "[$(date)] DONE Phase A"
ls -lh $OUT_DIR/
