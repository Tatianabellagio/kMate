#!/bin/bash
#SBATCH --job-name=haplo_v3qc_v3
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/haplo_v3qc_v3_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/haplo_v3qc_v3_%j.err

set -euo pipefail
BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping
SRC=$BASE/data/v3qc_v3/founders_231_v3qc_v3.vcf.gz
OUT=$BASE/data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz
BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
BGZIP=/home/tbellagio/miniforge3/envs/pang/bin/bgzip
TABIX=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/tabix

[ -s "$SRC" ] || { echo "ERROR: missing $SRC"; exit 1; }
[ ! -s "$OUT" ] || { echo "[$(date)] $OUT exists — exiting"; exit 0; }

echo "[$(date)] haploidize biallelic v3qc_v3 → $OUT"
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
