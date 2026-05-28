#!/bin/bash
#SBATCH --job-name=v3qc_merge
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/v3qc_merge_%j.out
#SBATCH --error=logs/v3qc_merge_%j.err

# =============================================================================
# build_v3qc_merged.sh
# v3qc Phase D: cactus_78 + pangenie_153_qc → founders_231_v3qc.vcf.gz
#   1. Subset cactus_pang69_1001g → cactus_78 (drop 5772 + 9947 samples)
#   2. Apply MAC<2 filter to PG side only (xwu's intent: drop short-read singletons
#      in inbred Arabidopsis). NOT applied to cactus side: cactus singletons are
#      genuine private-to-ecotype long-read variants, not calling artifacts.
#   3. Merge cactus_78 + pangenie_153_qc_mac2 on sample axis
#   4. fill-tags on merged panel
#   5. Drop AC=0 (catches ghost mono-REF records from removing 5772/9947 on cactus
#      side that also had no PG carriers — useless all-zero cn_var columns)
# Requires: pangenie_153_qc.vcf.gz from qc_pg_v4_filter.sh
# =============================================================================
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/hapfire_sv/pangenie_genotyping
MERGED=$BASE/data/merged
V3QC=$BASE/data/v3qc
TMP=$BASE/data/v3qc_tmp
mkdir -p $V3QC $TMP

BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/tabix

CACTUS80=$MERGED/cactus_pang69_1001g.vcf.gz
PG153_QC=$V3QC/pangenie_153_qc.vcf.gz

[ -s "$CACTUS80" ] || { echo "ERROR: missing $CACTUS80"; exit 1; }
[ -s "$PG153_QC" ] || { echo "ERROR: missing $PG153_QC — run qc_pg_v4_filter.sh first"; exit 1; }

# ---- Step 1: cactus_78 (drop 5772 + 9947) ----
echo "[$(date)] Step 1: cactus_80 → cactus_78 (drop 5772, 9947)"
CACTUS78=$TMP/cactus_78.vcf.gz
if [ ! -s "$CACTUS78" ]; then
    $BCF view -s ^5772,9947 $CACTUS80 --threads 8 -Oz -o $CACTUS78
    $TABIX -p vcf $CACTUS78
fi
N=$($BCF query -l $CACTUS78 | wc -l)
echo "  samples: $N (expected 78)"
[ "$N" = "78" ] || { echo "ERROR: got $N, expected 78"; exit 1; }

# ---- Step 1b: PG-side MAC<2 filter ----
# xwu's MAC>=2 is for inbred-Arabidopsis short-read SNP calls, where singletons
# are typically artifacts. Apply to PG side ONLY (long-read cactus singletons
# are real private alleles, kept). MAC<2 = AC<2 OR (AN-AC)<2 (symmetric).
echo "[$(date)] Step 1b: drop MAC<2 on PG side (153-sample subset) → pangenie_153_qc_mac2"
PG153_MAC2=$TMP/pangenie_153_qc_mac2.vcf.gz
if [ ! -s "$PG153_MAC2" ]; then
    N_PG_BEFORE=$($BCF index -n $PG153_QC)
    $BCF view -e 'INFO/AC<2 || (INFO/AN-INFO/AC)<2' $PG153_QC --threads 8 -Oz -o $PG153_MAC2
    $TABIX -p vcf $PG153_MAC2
    N_PG_AFTER=$($BCF index -n $PG153_MAC2)
    echo "  PG records before MAC<2: $N_PG_BEFORE"
    echo "  PG records after MAC<2:  $N_PG_AFTER"
    echo "  dropped: $((N_PG_BEFORE - N_PG_AFTER)) ($(awk -v b=$N_PG_BEFORE -v a=$N_PG_AFTER 'BEGIN{printf "%.2f%%", 100*(b-a)/b}'))"
fi

# ---- Step 2: merge cactus_78 + pangenie_153_qc_mac2 ----
echo "[$(date)] Step 2: merge cactus_78 + pangenie_153_qc_mac2 → founders_231_v3qc_raw"
MERGED_RAW=$TMP/founders_231_v3qc_raw.vcf.gz
if [ ! -s "$MERGED_RAW" ]; then
    $BCF merge $CACTUS78 $PG153_MAC2 --threads 8 -Oz -o $MERGED_RAW
    $TABIX -p vcf $MERGED_RAW
fi
N=$($BCF query -l $MERGED_RAW | wc -l)
echo "  samples: $N (expected 231)"
[ "$N" = "231" ] || { echo "ERROR: got $N, expected 231"; exit 1; }
N_REC_RAW=$($BCF index -n $MERGED_RAW)
echo "  records before MAC filter: $N_REC_RAW"

# ---- Step 3: fill-tags on merged panel ----
echo "[$(date)] Step 3: fill-tags on merged panel"
MERGED_FILLED=$TMP/founders_231_v3qc_filled.vcf.gz
if [ ! -s "$MERGED_FILLED" ]; then
    # NOTE: bcftools +fill-tags does NOT support MAC tag. Compute inline in filter.
    $BCF +fill-tags $MERGED_RAW --threads 8 -Oz -o $MERGED_FILLED -- -t AC,AN,AC_Het,F_MISSING
    $TABIX -p vcf $MERGED_FILLED
fi

# ---- Step 4: drop ONLY AC=0 records ("ghost" mono-REF from sample drop) ----
# REVISED 2026-05-16: previous MAC<2 filter was too aggressive — it dropped
# singletons (AC=1), which in a cactus-pangenome SV/SNP catalog are typically
# GENUINE private-to-ecotype alleles, not calling artifacts. xwu's MAC>=2
# applies to the GrENE-Net SNP-only context (inbred Arabidopsis short-read SNPs);
# our context (cactus SVs + multi-ploidy panel) needs singletons preserved.
#
# Singleton filtering (if needed for h-estimation) is better applied at the
# k-mer level in cn_full via per_sample_per_chrom.py --filt2 flag — NOT at
# the variant level in cn_var.
#
# Keep: any record with AC>=1 (at least one founder carries ALT)
# Drop: AC=0 records (mono-REF ghosts after sample drop — all-zero cn_var
#       columns, contribute nothing to EM)
echo "[$(date)] Step 4: drop AC=0 only (keep all singletons — they're real private-to-ecotype variants)"
FINAL=$V3QC/founders_231_v3qc.vcf.gz
$BCF view -e 'INFO/AC=0' $MERGED_FILLED --threads 8 -Oz -o $FINAL
$TABIX -p vcf $FINAL

N_REC_FINAL=$($BCF index -n $FINAL)
echo "  records after MAC<2 filter: $N_REC_FINAL"
echo "  dropped: $((N_REC_RAW - N_REC_FINAL)) ($(awk -v b=$N_REC_RAW -v a=$N_REC_FINAL 'BEGIN{printf "%.2f%%", 100*(b-a)/b}'))"

echo "[$(date)] DONE — founders_231_v3qc.vcf.gz built"
ls -lh $FINAL
