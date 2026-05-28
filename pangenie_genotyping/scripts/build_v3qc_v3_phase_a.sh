#!/bin/bash
#SBATCH --job-name=v3qc_v3_A
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=logs/v3qc_v3_A_%j.out
#SBATCH --error=logs/v3qc_v3_A_%j.err

# v3qc-v3 Phase A: replace V4 record-drop with per-cell het mask.
# Why: V4 dropped 27% of PG records (most with valid 1/1 carriers + one 0/1 noise call).
# Het mask preserves the hom-alt info, only NA's out the suspect 0/1 cells.
#
# ============================================================================
# ACTUAL DATA FLOW (read this carefully — there's a deceptive stale by-product)
# ============================================================================
#  pangenie_153_filled_bi.vcf.gz          (post-norm, pre-mask)
#       │
#       │ Step 1: +setGT het → ./.
#       ▼
#  pangenie_153_hetmasked_bi.vcf.gz       (intermediate)
#       │
#       │ Step 2: +fill-tags (recompute AC/AN/AC_Het/F_MISSING)
#       ▼
#  pangenie_153_hetmasked_filled_bi.vcf.gz   ←─── USED DOWNSTREAM (haploidize reads this)
#       │
#       │ Step 3: drop AC=0
#       ▼
#  pangenie_153_qc_v3.vcf.gz              ←─── DEAD END. NOT USED DOWNSTREAM.
#                                                 Kept only for historical/QC inspection.
#                                                 The per-side AC=0 drop loses real
#                                                 PG-REF info at cactus-private variants;
#                                                 that's why haploidize bypasses this file.
# ============================================================================
#
# Pipeline:
#   1. Reuse cactus_78_bi (already built)
#   2. Reuse pangenie_153_filled_bi (already built — post-norm + post-fill-tags)
#   3. bcftools +setGT to convert all 0/1 -> ./.  →  pangenie_153_hetmasked_bi
#   4. Re-fill-tags (AC, AN change after mask)    →  pangenie_153_hetmasked_filled_bi
#      ↑↑↑ THIS is the file that feeds the rest of the pipeline (via haploidize) ↑↑↑
#   5. Drop AC=0 records                          →  pangenie_153_qc_v3 (DEAD-END, NOT USED)
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/hapfire_sv/pangenie_genotyping/data
BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/tabix

CACTUS78_BI=$BASE/v3qc_v2/cactus_78_bi.vcf.gz       # reused
PG_FILLED=$BASE/v3qc_v2/pangenie_153_filled_bi.vcf.gz  # reused

OUT_DIR=$BASE/v3qc_v3
mkdir -p $OUT_DIR

PG_HETMASKED=$OUT_DIR/pangenie_153_hetmasked_bi.vcf.gz
PG_HETMASKED_FILLED=$OUT_DIR/pangenie_153_hetmasked_filled_bi.vcf.gz
PG_QC=$OUT_DIR/pangenie_153_qc_v3.vcf.gz

[ -s "$CACTUS78_BI" ] || { echo "ERROR: missing $CACTUS78_BI"; exit 1; }
[ -s "$PG_FILLED" ] || { echo "ERROR: missing $PG_FILLED"; exit 1; }

# Reuse cactus_78_bi by symlink
ln -sf $CACTUS78_BI $OUT_DIR/cactus_78_bi.vcf.gz
ln -sf $CACTUS78_BI.tbi $OUT_DIR/cactus_78_bi.vcf.gz.tbi

# Step 1: het mask (convert all 0/1 -> ./.)
echo "[$(date)] Step 1: het mask (set 0/1 -> ./.)"
N_BEFORE=$($BCF index -n $PG_FILLED)
echo "  PG records pre-mask: $N_BEFORE"
if [ ! -s "$PG_HETMASKED" ]; then
    $BCF +setGT $PG_FILLED --threads 8 -Oz -o $PG_HETMASKED -- -t q -i 'GT="het"' -n .
    $TABIX -p vcf $PG_HETMASKED
fi

# Step 2: re-fill-tags (AC, AN, AC_Het, F_MISSING all change after masking)
echo "[$(date)] Step 2: re-fill-tags AC, AN, AC_Het, F_MISSING"
if [ ! -s "$PG_HETMASKED_FILLED" ]; then
    $BCF +fill-tags $PG_HETMASKED --threads 8 -Oz -o $PG_HETMASKED_FILLED -- -t AC,AN,AC_Het,F_MISSING
    $TABIX -p vcf $PG_HETMASKED_FILLED
fi

# Verify het mask worked: AC_Het should be 0 everywhere now
N_HET_REMAINING=$($BCF view -H -e 'INFO/AC_Het=0' $PG_HETMASKED_FILLED 2>/dev/null | wc -l)
echo "  Records with AC_Het>0 after masking (should be 0): $N_HET_REMAINING"

# Step 3: drop records with AC=0 after mask.
#
# WARNING: This output ($PG_QC = pangenie_153_qc_v3.vcf.gz) is a DEAD-END.
# It is NOT used by haploidize_pg_hetmasked.sh or by Phase B. The haploidize
# step reads pangenie_153_hetmasked_filled_bi.vcf.gz directly (the PRE-AC=0
# file), so this AC=0 cleanup never affects the production panel.
#
# Why we don't apply the AC=0 drop to the production chain:
#   A record where PG has all-REF (0/0) cells but no PG carriers gets AC=0
#   after het mask. Dropping it would LOSE the real PG REF information at
#   cactus-private variants. The proper place for the AC=0 cleanup is
#   POST-MERGE (Phase B Step 4), where AC=0 means "no carriers anywhere
#   across the whole 231 panel."
#
# The Phase A AC=0 output is still produced for historical/QC inspection
# (e.g. counting how many records would have been lost at various stages).
# Don't use it as input anywhere downstream.
echo "[$(date)] Step 3: drop AC=0 → pangenie_153_qc_v3 (DEAD-END, NOT USED DOWNSTREAM)"
$BCF view -e 'INFO/AC=0' $PG_HETMASKED_FILLED --threads 8 -Oz -o $PG_QC
$TABIX -p vcf $PG_QC
N_AFTER=$($BCF index -n $PG_QC)
DROPPED=$((N_BEFORE - N_AFTER))
PCT=$(awk -v b=$N_BEFORE -v a=$N_AFTER 'BEGIN{printf "%.2f%%", 100*(b-a)/b}')
echo "  PG records after het-mask + AC=0 drop:  $N_AFTER  (NOT USED DOWNSTREAM)"
echo "  Total dropped: $DROPPED ($PCT)"
echo "  NOTE: Downstream pipeline uses pangenie_153_hetmasked_filled_bi (pre-AC=0)."

echo ""
echo "[$(date)] DONE Phase A"
ls -lh $OUT_DIR/
