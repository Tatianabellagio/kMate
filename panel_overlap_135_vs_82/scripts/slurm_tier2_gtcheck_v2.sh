#!/bin/bash
#SBATCH --job-name=overlap_t2_gtchk2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/tier2_gtcheck_v2_%j.out
#SBATCH --error=logs/tier2_gtcheck_v2_%j.err
mkdir -p logs
set -euo pipefail

# Tier 2 v2: convert haploid cactus GT to diploid (0 → 0/0, 1 → 1/1, . → ./.),
# write a new VCF, then run bcftools gtcheck.

source /global/home/users/tbellg/miniforge3/etc/profile.d/conda.sh
conda activate sequencing_pipeline

WORK=/global/scratch/users/tbellg/hapfire_sv/panel_overlap_135_vs_82
SNPS=$WORK/data/pang135_biallelic_snps.vcf.gz
DIP=$WORK/data/pang135_biallelic_snps_diploid.vcf.gz

if [ ! -s "$DIP" ]; then
  echo "[$(date +%H:%M:%S)] converting haploid GT to diploid"
  # zcat header lines verbatim; for data lines, transform fields 10+ (GT) by mapping 0→0/0, 1→1/1, .→./.
  zcat "$SNPS" \
    | awk -v OFS='\t' '
        /^#/ {print; next}
        {
          for (i=10; i<=NF; i++) {
            if ($i == "0")      $i = "0/0";
            else if ($i == "1") $i = "1/1";
            else if ($i == ".") $i = "./.";
            # else: leave as-is (already diploid?)
          }
          print
        }' \
    | bgzip -@ 4 > "$DIP"
  tabix -p vcf -f "$DIP"
  echo "[$(date +%H:%M:%S)] diploid SNP file: $(stat -c%s $DIP) bytes"
fi

PAIRS=$WORK/data/gtcheck_pairs.tsv
OUT=$WORK/results/tier2_gtcheck.tsv

echo "[$(date +%H:%M:%S)] running gtcheck on $(wc -l < $PAIRS) pairs"
bcftools gtcheck --use GT --no-HWE-prob -P $PAIRS -O t -o $OUT $DIP
echo "[$(date +%H:%M:%S)] gtcheck done → $OUT"
echo "[$(date +%H:%M:%S)] head:"
head -25 $OUT
