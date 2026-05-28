#!/bin/bash
#SBATCH --job-name=haplo_pg_hm
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/haplo_pg_hm_%j.out
#SBATCH --error=logs/haplo_pg_hm_%j.err

# Pre-haploidize pangenie_153_hetmasked_filled_bi.vcf.gz BEFORE the cactus merge.
# Reason: bcftools merge of diploid PG + haploid cactus_78 has a context-dependent
# bug that drops PG GTs at AC=0 records (manifests on full-dataset merges of
# ~7.7M records; not on smaller subsets). Pre-haploidizing PG eliminates the
# ploidy mismatch, so merge sees haploid+haploid uniformly.
#
# Diploid → haploid: 0/0 → 0, 1/1 → 1, ./. → ., else → . (paranoid fallback).
mkdir -p logs
set -euo pipefail
BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping
SRC=$BASE/data/v3qc_v3/pangenie_153_hetmasked_filled_bi.vcf.gz
OUT=$BASE/data/v3qc_v3/pangenie_153_hetmasked_haploid.vcf.gz
BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
BGZIP=/global/home/users/tbellg/miniforge3/envs/pang/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/tabix

[ -s "$SRC" ] || { echo "ERROR: missing $SRC"; exit 1; }
[ ! -s "$OUT" ] || { echo "[$(date)] $OUT exists — exiting"; exit 0; }

echo "[$(date)] haploidize PG het-masked → $OUT"
$BCF view --threads 2 -Ov $SRC | \
awk 'BEGIN{OFS="\t"}
    /^##/ { print; next }
    /^#CHROM/ { print; next }
    {
        # Keep INFO (col 8) but rebuild FORMAT (col 9) and sample fields.
        # Only need GT in the output — drop GQ/GL/KC for compactness.
        $9 = "GT"
        for (i=10; i<=NF; i++) {
            split($i, parts, ":")
            g = parts[1]
            if (g == "0/0" || g == "0|0") { $i = "0" }
            else if (g == "1/1" || g == "1|1") { $i = "1" }
            else if (g == "./." || g == ".|." || g == ".") { $i = "." }
            else { $i = "." }    # paranoid fallback for any unexpected GT
        }
        print
    }' | $BGZIP -@ 4 -c > $OUT
$TABIX -p vcf $OUT

N_REC=$($BCF index -n $OUT)
N_SAM=$($BCF query -l $OUT | wc -l)
echo "[$(date)] DONE records=$N_REC samples=$N_SAM"
ls -lh $OUT
