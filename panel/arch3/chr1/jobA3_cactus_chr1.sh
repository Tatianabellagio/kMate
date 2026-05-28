#!/bin/bash
#SBATCH --job-name=chr1_cactus
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=logs/A3_cactus_%j.out
#SBATCH --error=logs/A3_cactus_%j.err
mkdir -p logs
set -euo pipefail

# Phase 2 Job A3: cactus-side full pipeline on full Chr1 (depends on A1).
# Steps: subset cactus_78 → transfer_id → convert-to-biallelic → sort/bgzip (already haploid).

cd /global/scratch/users/tbellg/kmate/panel/arch3/chr1

BCF=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bcftools
BGZIP=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/gwas/bin/tabix
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CONVERT=/global/scratch/users/tbellg/kmate/external_tools/pangenie-tools/pipelines/run-from-callset/scripts/convert-to-biallelic.py
TRANSFER=/global/scratch/users/tbellg/kmate/scratch/arch3_test/transfer_id_annotation.py

CHR1_ANNOT=chr1_135_annotated.sorted.vcf.gz
BIAL_CATALOG=chr1_135_annotated_biallelic.sorted.vcf.gz
CACTUS78_RAW=/global/scratch/users/tbellg/kmate/panel/pangenie_genotyping/data/v3qc_tmp/cactus_78.vcf.gz

[ -s $CHR1_ANNOT ] || { echo "ERROR: missing $CHR1_ANNOT (A1 not done)"; exit 1; }

echo "[$(date)] === Step 1: subset cactus_78 to full Chr1 ==="
CACTUS_CHR1=cactus_78_chr1.vcf.gz
if [ ! -s $CACTUS_CHR1 ]; then
  $BCF view -r Chr1 $CACTUS78_RAW -Oz -o $CACTUS_CHR1
  $TABIX -p vcf $CACTUS_CHR1
fi
echo "  cactus_78 Chr1 records: $($BCF view -H $CACTUS_CHR1 | wc -l)"
echo "  samples: $($BCF query -l $CACTUS_CHR1 | wc -l)"

echo
echo "[$(date)] === Step 2: transfer INFO/ID from chr1 catalog ==="
CACTUS_ANNOT=cactus_78_chr1_annotated.vcf
$PY $TRANSFER --cactus $CHR1_ANNOT --pg $CACTUS_CHR1 --out $CACTUS_ANNOT

echo
echo "[$(date)] === Step 3: convert-to-biallelic ==="
CACTUS_BIAL=cactus_78_chr1_biallelic.vcf
cat $CACTUS_ANNOT | $PY $CONVERT $BIAL_CATALOG > $CACTUS_BIAL 2> cactus_convert.log
echo "  exit: $?"
echo "  records emitted: $(grep -vc '^#' $CACTUS_BIAL)"

echo
echo "[$(date)] === Step 4: sort + bgzip + index cactus_78 biallelic (already haploid) ==="
CACTUS_HAP=cactus_78_chr1_haploid.vcf.gz
awk '$1 ~ /^#/ {print $0; next} {print $0 | "sort -k1,1 -k2,2n"}' $CACTUS_BIAL > cactus_78_chr1_biallelic.sorted.vcf
$BGZIP -f -c cactus_78_chr1_biallelic.sorted.vcf > $CACTUS_HAP
$TABIX -p vcf $CACTUS_HAP
echo "  $CACTUS_HAP: $(ls -lh $CACTUS_HAP | awk '{print $5}')"

echo
echo "[$(date)] DONE A3 (cactus Chr1)"
