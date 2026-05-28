#!/bin/bash
#SBATCH --job-name=haplo_v3qc_v2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/haplo_v3qc_v2_%j.out
#SBATCH --error=logs/haplo_v3qc_v2_%j.err

# Haploidize founders_231_v3qc_v2.vcf.gz → founders_231_v3qc_v2.haploid.vcf.gz
# Already biallelic, so we just need the AWK haploidization (no norm step needed).
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/kmate/panel/pangenie_genotyping
SRC=$BASE/data/v3qc_v2/founders_231_v3qc_v2.vcf.gz
OUT=$BASE/data/v3qc_v2/founders_231_v3qc_v2.haploid.vcf.gz
BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
BGZIP=/global/home/users/tbellg/miniforge3/envs/pang/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/tabix

[ -s "$SRC" ] || { echo "ERROR: missing $SRC"; exit 1; }
[ ! -s "$OUT" ] || { echo "[$(date)] $OUT exists — exiting"; exit 0; }

echo "[$(date)] haploidize biallelic v3qc_v2 → $OUT"
$BCF view --threads 2 -Ov $SRC | \
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
    }' | $BGZIP -@ 4 -c > $OUT
$TABIX -p vcf $OUT

N_REC=$($BCF index -n $OUT)
N_SAM=$($BCF query -l $OUT | wc -l)
echo "[$(date)] DONE records=$N_REC samples=$N_SAM"
ls -lh $OUT
