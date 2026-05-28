#!/bin/bash
#SBATCH --job-name=lowmiss_vcf
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/lowmiss_vcf_%j.out
#SBATCH --error=logs/lowmiss_vcf_%j.err

# Lenient version of the no-missing test: drop records with F_MISSING > 0.5.
# Catches the PG-MAC merge artifact (F_MISSING = 0.66) but keeps records with a few GQ-masked cells.
# Expected: keeps ~32% of records (vs 18% for strict F_MISSING==0).
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/kmate/panel/pangenie_genotyping
SRC=$BASE/data/v3qc/founders_231_v3qc.haploid.vcf.gz
OUT=$BASE/data/v3qc/founders_231_v3qc.lowmiss.haploid.vcf.gz

BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/tabix

[ -s "$SRC" ] || { echo "ERROR: missing $SRC"; exit 1; }
[ ! -s "$OUT" ] || { echo "[$(date)] $OUT exists — exiting"; exit 0; }

echo "[$(date)] Recompute F_MISSING and drop records with F_MISSING > 0.5"
$BCF +fill-tags $SRC --threads 4 -- -t F_MISSING 2>/dev/null | \
    $BCF view -e 'INFO/F_MISSING > 0.5' --threads 4 -Oz -o $OUT
$TABIX -p vcf $OUT

N_BEFORE=$($BCF index -n $SRC)
N_AFTER=$($BCF index -n $OUT)
echo "[$(date)] DONE"
echo "  records before: $N_BEFORE"
echo "  records after:  $N_AFTER"
echo "  dropped: $((N_BEFORE - N_AFTER)) ($(awk -v b=$N_BEFORE -v a=$N_AFTER 'BEGIN{printf "%.2f%%", 100*(b-a)/b}'))"
ls -lh $OUT
