#!/bin/bash
#SBATCH --job-name=haplo_v3qc
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=logs/haplo_v3qc_%j.out
#SBATCH --error=logs/haplo_v3qc_%j.err

# =============================================================================
# haploidize_v3qc_vcf.sh
# Build founders_231_v3qc.haploid.vcf.gz from founders_231_v3qc.vcf.gz.
# Carrier-status haploid: 0/0 → 0, ./. → ., any-ALT → 1.
# Mirrors haploidize_merged_vcf.sh for the v3qc input.
# =============================================================================
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping
SRC=$BASE/data/v3qc/founders_231_v3qc.vcf.gz
OUT=$BASE/data/v3qc/founders_231_v3qc.haploid.vcf.gz

BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/tabix

[ -s "$SRC" ] || { echo "ERROR: missing $SRC" >&2; exit 1; }
[ ! -s "$OUT" ] || { echo "[$(date)] $OUT exists — exiting (delete to rerun)"; exit 0; }

echo "[$(date)] decompose multi-allelics → haploidize → bgzip"

$BCF norm -m -any --threads 4 -Ou "$SRC" 2>/dev/null | \
$BCF view --threads 2 -Ov - | \
awk 'BEGIN{OFS="\t"}
    /^##/ { print; next }
    /^#CHROM/ { print; next }
    {
        $9 = "GT"
        for (i=10; i<=NF; i++) {
            split($i, parts, ":")
            g = parts[1]
            if (g == "0" || g == "0/0" || g == "0|0") {
                $i = "0"
            } else if (g == "." || g == "./." || g == ".|.") {
                $i = "."
            } else {
                $i = "1"
            }
        }
        print
    }' | \
/global/home/users/tbellg/miniforge3/envs/pang/bin/bgzip -@ 4 -c > "$OUT"

$TABIX -p vcf "$OUT"

N_REC=$($BCF view -H "$OUT" 2>/dev/null | wc -l)
N_SAM=$($BCF query -l "$OUT" 2>/dev/null | wc -l)
echo "[$(date)] DONE  records=$N_REC  samples=$N_SAM"
ls -lh "$OUT"
