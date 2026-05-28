#!/bin/bash
#SBATCH --job-name=nomac_f50
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --output=logs/nomac_f50_%j.out
#SBATCH --error=logs/nomac_f50_%j.err

# Build a v3qc-style merged + haploid VCF SKIPPING PG-MAC step and applying F_MISSING<=0.5.
# Pipeline:
#   merge cactus_78 + pangenie_153_qc (no PG-MAC filter)
#   fill-tags AC, AN, F_MISSING
#   drop records where AC=0 OR F_MISSING > 0.5
#   norm -m -any + haploidize (same convention as v3 / v3qc)
# Output: founders_231_v3qc_noMAC_F50.haploid.vcf.gz
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping
BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
BGZIP=/global/home/users/tbellg/miniforge3/envs/pang/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/tabix

CACTUS78=$BASE/data/v3qc_tmp/cactus_78.vcf.gz
PG_QC=$BASE/data/v3qc/pangenie_153_qc.vcf.gz
TMP=$BASE/data/v3qc_tmp
OUT_RAW=$TMP/founders_231_v3qc_noMAC_F50.vcf.gz
OUT_HAP=$BASE/data/v3qc/founders_231_v3qc_noMAC_F50.haploid.vcf.gz

[ -s "$CACTUS78" ] || { echo "ERROR: missing $CACTUS78"; exit 1; }
[ -s "$PG_QC" ] || { echo "ERROR: missing $PG_QC"; exit 1; }
[ ! -s "$OUT_HAP" ] || { echo "[$(date)] $OUT_HAP exists — exiting"; exit 0; }

# Step 1: merge + fill-tags + drop AC=0 OR F_MISSING>0.5
echo "[$(date)] Step 1: merge cactus_78 + pangenie_153_qc, fill-tags, drop AC=0 OR F_MISSING>0.5"
if [ ! -s "$OUT_RAW" ]; then
    $BCF merge $CACTUS78 $PG_QC --threads 8 -Ou | \
        $BCF +fill-tags --threads 8 -- -t AC,AN,F_MISSING | \
        $BCF view -e 'INFO/AC=0 || INFO/F_MISSING > 0.5' --threads 8 -Oz -o $OUT_RAW
    $TABIX -p vcf $OUT_RAW
fi
N=$($BCF query -l $OUT_RAW | wc -l)
N_REC=$($BCF index -n $OUT_RAW)
echo "  samples: $N (expected 231)"
echo "  records: $N_REC"
[ "$N" = "231" ] || { echo "ERROR: got $N samples"; exit 1; }

# Step 2: norm -m -any + haploidize (same as haploidize_v3qc_vcf.sh)
echo "[$(date)] Step 2: norm + haploidize → $OUT_HAP"
$BCF norm -m -any --threads 4 -Ou "$OUT_RAW" 2>/dev/null | \
$BCF view --threads 2 -Ov - | \
awk 'BEGIN{OFS="\t"}
    /^##/ { print; next }
    /^#CHROM/ { print; next }
    {
        $9 = "GT"
        for (i=10; i<=NF; i++) {
            split($i, parts, ":")
            g = parts[1]
            if (g == "0" || g == "0/0" || g == "0|0") { $i = "0" }
            else if (g == "." || g == "./." || g == ".|.") { $i = "." }
            else { $i = "1" }
        }
        print
    }' | $BGZIP -@ 4 -c > $OUT_HAP
$TABIX -p vcf $OUT_HAP

N_HAP=$($BCF index -n $OUT_HAP)
echo "[$(date)] DONE"
echo "  haploid records: $N_HAP"
ls -lh $OUT_HAP
