#!/bin/bash
#SBATCH --job-name=v3qc_v3_B
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/v3qc_v3_B_%j.out
#SBATCH --error=logs/v3qc_v3_B_%j.err

# v3qc-v3 Phase B: merge cactus_78_bi (haploid) + pangenie_153_hetmasked_haploid → re-decompose → AC=0 cleanup.
# Same logic as v3qc_v2 Phase B but on the het-masked PG side.
#
# ============================================================================
# IMPORTANT: PG input is pangenie_153_hetmasked_haploid.vcf.gz (NOT pangenie_153_qc_v3).
# pangenie_153_qc_v3.vcf.gz is a dead-end Phase A by-product with the AC=0 drop applied;
# we DON'T want that AC=0 drop on the PG side because it loses real PG-REF info at
# cactus-private variants. The AC=0 cleanup is properly done POST-MERGE (Step 4 below).
# See header comment in build_v3qc_v3_phase_a.sh for full data-flow diagram.
# ============================================================================
mkdir -p logs
set -euo pipefail

# Repo dir for this stage; override $PANGENIE_GT for sbatch spool copies.
BASE="${PANGENIE_GT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"/data
BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/tabix
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa

CACTUS78_BI=$BASE/v3qc_v3/cactus_78_bi.vcf.gz
# PG-side input: the haploidized het-masked PG file.
# - "het-masked" : 0/1 cells → ./. (per Phase A step 1)
# - "haploid"    : ./.→. , 0/0→0, 1/1→1 (via haploidize_pg_hetmasked.sh awk)
# - NOT post-AC=0-filtered (so PG REF info at cactus-private variants is preserved).
# Pre-haploidizing PG before merge also avoids a bcftools merge bug with
# mixed-ploidy + AC=0 + large record count that drops PG GTs at AC=0 records.
# Verified empirically on 11 test positions (Chr1:751..1898).
PG_HETMASKED_HAPLOID=$BASE/v3qc_v3/pangenie_153_hetmasked_haploid.vcf.gz
TMP_DIR=$BASE/v3qc_v3/tmp
mkdir -p $TMP_DIR
OUT_DIR=$BASE/v3qc_v3
MERGED_RAW=$TMP_DIR/founders_231_v3qc_v3_merged_raw.vcf.gz
MERGED_BI=$TMP_DIR/founders_231_v3qc_v3_merged_bi.vcf.gz
MERGED_FILLED=$TMP_DIR/founders_231_v3qc_v3_merged_filled.vcf.gz
FINAL=$OUT_DIR/founders_231_v3qc_v3.vcf.gz

[ -s "$CACTUS78_BI" ] || { echo "ERROR: missing $CACTUS78_BI"; exit 1; }
[ -s "$PG_HETMASKED_HAPLOID" ] || { echo "ERROR: missing $PG_HETMASKED_HAPLOID"; exit 1; }
[ ! -s "$FINAL" ]    || { echo "[$(date)] $FINAL exists — exiting"; exit 0; }

# Step 1: merge
echo "[$(date)] Step 1: merge cactus_78_bi + pangenie_153_hetmasked_haploid → merged_raw"
if [ ! -s "$MERGED_RAW" ]; then
    $BCF merge $CACTUS78_BI $PG_HETMASKED_HAPLOID --threads 8 -Oz -o $MERGED_RAW
    $TABIX -p vcf $MERGED_RAW
fi
N_SAM=$($BCF query -l $MERGED_RAW | wc -l)
N_REC=$($BCF index -n $MERGED_RAW)
echo "  samples: $N_SAM (expected 231)"
echo "  records: $N_REC"
[ "$N_SAM" = "231" ] || { echo "ERROR: got $N_SAM samples"; exit 1; }

# Step 2: re-decompose
echo "[$(date)] Step 2: re-decompose merged VCF"
if [ ! -s "$MERGED_BI" ]; then
    $BCF norm -f $REF -m -any $MERGED_RAW --threads 8 -Ou 2> $TMP_DIR/merged_norm.log | \
        $BCF sort -m 8G -Oz -o $MERGED_BI -
    $TABIX -p vcf $MERGED_BI
fi
N_BI=$($BCF index -n $MERGED_BI)
echo "  records after re-decompose: $N_BI"

# Step 3: fill-tags
echo "[$(date)] Step 3: fill-tags"
if [ ! -s "$MERGED_FILLED" ]; then
    $BCF +fill-tags $MERGED_BI --threads 8 -Oz -o $MERGED_FILLED -- -t AC,AN,AC_Het,F_MISSING
    $TABIX -p vcf $MERGED_FILLED
fi

# Step 4: AC=0 cleanup (POST-MERGE — this is the ONLY place AC=0 should be applied)
# At this point AC=0 means "no carriers in the WHOLE 231 panel" — truly variantless.
# Safe to drop because the record has no useful AF signal for downstream EM.
echo "[$(date)] Step 4: drop AC=0 (post-merge: AC=0 means variantless across 231 panel)"
N_PRE=$($BCF index -n $MERGED_FILLED)
$BCF view -e 'INFO/AC=0' $MERGED_FILLED --threads 8 -Oz -o $FINAL
$TABIX -p vcf $FINAL
N_FINAL=$($BCF index -n $FINAL)
echo "  records before AC=0: $N_PRE"
echo "  records after AC=0:  $N_FINAL"

echo ""
echo "=== F_MISSING distribution on founders_231_v3qc_v3 (Chr1) ==="
$BCF query -r Chr1 -f '%INFO/F_MISSING\n' $FINAL | \
/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -c "
import sys, numpy as np
fm = np.array([float(x.strip()) for x in sys.stdin if x.strip()])
N = len(fm)
print(f'Total Chr1 records: {N:,}')
buckets = [
    ('F_MISSING == 0',                fm == 0),
    ('0 < F <= 0.05',                 (fm > 0) & (fm <= 0.05)),
    ('0.05 < F <= 0.20',              (fm > 0.05) & (fm <= 0.20)),
    ('0.20 < F <= 0.50',              (fm > 0.20) & (fm <= 0.50)),
    ('0.50 < F <= 0.80',              (fm > 0.50) & (fm <= 0.80)),
    ('F > 0.80',                      fm > 0.80),
]
for label, mask in buckets:
    n = mask.sum()
    print(f'  {label:30s}: {n:>10,} ({100*n/N:>5.2f}%)')
"

echo ""
echo "[$(date)] DONE Phase B"
ls -lh $FINAL
