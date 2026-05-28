#!/bin/bash
#SBATCH --job-name=v3qc_v2_B
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/v3qc_v2_B_%j.out
#SBATCH --error=logs/v3qc_v2_B_%j.err

# Phase B — merge biallelic cactus_78 + pangenie_153_qc_v2 → re-decompose → AC=0 cleanup.
# Also reports F_MISSING distribution on the final merged VCF for evaluation.
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/hapfire_sv/pangenie_genotyping/data
BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/tabix
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa

CACTUS78_BI=$BASE/v3qc_v2/cactus_78_bi.vcf.gz
PG_QC_V2=$BASE/v3qc_v2/pangenie_153_qc_v2.vcf.gz
TMP_DIR=$BASE/v3qc_v2/tmp
mkdir -p $TMP_DIR
OUT_DIR=$BASE/v3qc_v2
MERGED_RAW=$TMP_DIR/founders_231_v3qc_v2_merged_raw.vcf.gz
MERGED_BI=$TMP_DIR/founders_231_v3qc_v2_merged_bi.vcf.gz
MERGED_FILLED=$TMP_DIR/founders_231_v3qc_v2_merged_filled.vcf.gz
FINAL=$OUT_DIR/founders_231_v3qc_v2.vcf.gz

[ -s "$CACTUS78_BI" ] || { echo "ERROR: missing $CACTUS78_BI"; exit 1; }
[ -s "$PG_QC_V2"   ] || { echo "ERROR: missing $PG_QC_V2"; exit 1; }
[ ! -s "$FINAL" ] || { echo "[$(date)] $FINAL exists — exiting"; exit 0; }

# Step 1: merge
echo "[$(date)] Step 1: merge cactus_78_bi + pangenie_153_qc_v2 → merged_raw"
if [ ! -s "$MERGED_RAW" ]; then
    $BCF merge $CACTUS78_BI $PG_QC_V2 --threads 8 -Oz -o $MERGED_RAW
    $TABIX -p vcf $MERGED_RAW
fi
N_SAM=$($BCF query -l $MERGED_RAW | wc -l)
N_REC=$($BCF index -n $MERGED_RAW)
echo "  samples: $N_SAM (expected 231)"
echo "  records: $N_REC"
[ "$N_SAM" = "231" ] || { echo "ERROR: got $N_SAM samples"; exit 1; }

# Step 2: re-decompose (catch multi-allelics created by merge at colocated positions)
echo "[$(date)] Step 2: re-decompose merged VCF (catch multi-allelics from merge)"
if [ ! -s "$MERGED_BI" ]; then
    $BCF norm -f $REF -m -any $MERGED_RAW --threads 8 -Ou 2> $TMP_DIR/merged_norm.log | \
        $BCF sort -m 8G -Oz -o $MERGED_BI -
    $TABIX -p vcf $MERGED_BI
fi
N_BI=$($BCF index -n $MERGED_BI)
echo "  records after re-decompose: $N_BI"

# Step 3: fill-tags on biallelic merged
echo "[$(date)] Step 3: fill-tags AC, AN, AC_Het, F_MISSING on biallelic merged"
if [ ! -s "$MERGED_FILLED" ]; then
    $BCF +fill-tags $MERGED_BI --threads 8 -Oz -o $MERGED_FILLED -- -t AC,AN,AC_Het,F_MISSING
    $TABIX -p vcf $MERGED_FILLED
fi

# Step 4: AC=0 cleanup (now per-ALT correct)
echo "[$(date)] Step 4: drop AC=0 (per-ALT correct on biallelic)"
N_PRE=$($BCF index -n $MERGED_FILLED)
$BCF view -e 'INFO/AC=0' $MERGED_FILLED --threads 8 -Oz -o $FINAL
$TABIX -p vcf $FINAL
N_FINAL=$($BCF index -n $FINAL)
echo "  records before AC=0: $N_PRE"
echo "  records after AC=0:  $N_FINAL"
echo "  dropped: $((N_PRE - N_FINAL)) ($(awk -v b=$N_PRE -v a=$N_FINAL 'BEGIN{printf "%.2f%%", 100*(b-a)/b}'))"

echo ""
echo "[$(date)] DONE Phase B"
ls -lh $FINAL

# Step 5: F_MISSING distribution (Chr1 only for speed; extrapolate to genome)
echo ""
echo "=== F_MISSING distribution on founders_231_v3qc_v2 (Chr1) ==="
$BCF query -r Chr1 -f '%INFO/F_MISSING\n' $FINAL | \
/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python -c "
import sys, numpy as np
fm = np.array([float(x.strip()) for x in sys.stdin if x.strip()])
N = len(fm)
print(f'Total Chr1 records: {N:,}')
buckets = [
    ('F_MISSING == 0 (no missings)',  fm == 0),
    ('0 < F <= 0.05',                 (fm > 0) & (fm <= 0.05)),
    ('0.05 < F <= 0.20',              (fm > 0.05) & (fm <= 0.20)),
    ('0.20 < F <= 0.50',              (fm > 0.20) & (fm <= 0.50)),
    ('0.50 < F <= 0.80',              (fm > 0.50) & (fm <= 0.80)),
    ('F > 0.80',                      fm > 0.80),
]
for label, mask in buckets:
    n = mask.sum()
    print(f'  {label:42s}: {n:>10,} ({100*n/N:>5.2f}%)')
"
