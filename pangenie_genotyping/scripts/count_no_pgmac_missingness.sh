#!/bin/bash
#SBATCH --job-name=count_nomac
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:30:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/count_nomac_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/count_nomac_%j.err

# Compute F_MISSING distribution if we drop PG-MAC<2:
#   merge cactus_78 (5.2M records) + pangenie_153_qc (3.5M, post-V4 only) → ~5.2M merged
#   +fill-tags F_MISSING on merged
#   Count records by F_MISSING bucket (Chr1 only for speed)
set -euo pipefail

BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/data

CACTUS78=$BASE/v3qc_tmp/cactus_78.vcf.gz
PG_QC=$BASE/v3qc/pangenie_153_qc.vcf.gz
OUT_TSV=$BASE/v3qc/no_pgmac_missingness_chr1.tsv

[ -s "$CACTUS78" ] || { echo "ERROR: missing $CACTUS78"; exit 1; }
[ -s "$PG_QC" ] || { echo "ERROR: missing $PG_QC"; exit 1; }

echo "[$(date)] Stream merge + fill-tags + extract F_MISSING for Chr1 only"

$BCF merge $CACTUS78 $PG_QC -r Chr1 --threads 4 2>/dev/null | \
  $BCF +fill-tags - --threads 4 -- -t F_MISSING,AC 2>/dev/null | \
  $BCF query -f '%INFO/F_MISSING\t%INFO/AC\n' > $OUT_TSV

echo "[$(date)] DONE: $(wc -l < $OUT_TSV) records → $OUT_TSV"

# Summarize
/home/tbellagio/miniforge3/envs/hapfm/bin/python << EOF
import numpy as np
data = np.loadtxt("$OUT_TSV", dtype=str)
fm = data[:,0].astype(float)
ac = data[:,1].astype(int)
N = len(fm)

print(f"\n=== F_MISSING distribution: merged cactus_78 + pangenie_153_qc (NO PG-MAC) — Chr1 ===")
print(f"Total records: {N:,}")
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

# AC=0 records in this merged-no-MAC: should be much fewer than the current v3qc state
ac0 = (ac == 0).sum()
print(f"\n  Records with AC=0 in merged: {ac0:,} ({100*ac0/N:.2f}%)")
print(f"  (these get dropped later by AC=0 cleanup)")

# Strict / lenient thresholds counts
strict = (fm == 0).sum()
lenient = (fm <= 0.5).sum()
print(f"\n=== If we apply F_MISSING filter on this merged VCF ===")
print(f"  Strict (F==0):  keep {strict:,} ({100*strict/N:.1f}%)")
print(f"  Lenient (F<=0.5): keep {lenient:,} ({100*lenient/N:.1f}%)")
print(f"\nExtrapolation to all 5 chroms (assume Chr1 ~25% of records):")
print(f"  Total merged records genome-wide: ~{N*4:,}")
print(f"  Strict survivors:  ~{strict*4:,}")
print(f"  Lenient survivors: ~{lenient*4:,}")
EOF

echo "[$(date)] DONE"
