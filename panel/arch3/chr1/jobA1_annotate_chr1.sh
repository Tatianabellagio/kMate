#!/bin/bash
#SBATCH --job-name=chr1_annot
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=logs/A1_annot_%j.out
#SBATCH --error=logs/A1_annot_%j.err
mkdir -p logs
set -euo pipefail

# Phase 2 Job A1: run annotate_vcf on FULL Chr1 (135-sample input).
# Produces the annotated multi-allelic VCF + biallelic catalog for Chr1.
# These outputs feed both PG (A2) and cactus (A3) downstream jobs.

# Run in this script's directory; override $ARCH3_CHR1_DIR when launching from
# an sbatch spool copy outside the source tree.
cd "${ARCH3_CHR1_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

BCF=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bcftools
BGZIP=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/gwas/bin/tabix
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
ANNOTATE=../../../external/genotyping-pipelines/prepare-vcf-MC/workflow/scripts/annotate_vcf.py  # see README Prerequisites (external/ is gitignored)

GFA=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.gfa.gz
FULL_135_VCF=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.vcf.gz

REGION="Chr1"   # full Chr1
OUTPREFIX=chr1_135_annotated

echo "[$(date)] === Step 1: subset 135-sample VCF to full Chr1 ==="
VCF_CHR1=chr1_135.vcf
if [ ! -s ${VCF_CHR1}.done ]; then
  $BCF view -r $REGION $FULL_135_VCF > $VCF_CHR1
  touch ${VCF_CHR1}.done
fi
N=$($BCF view -H $VCF_CHR1 | wc -l)
NS=$($BCF view -h $VCF_CHR1 | grep "^#CHROM" | awk '{print NF-9}')
echo "  records: $N (full Chr1)"
echo "  samples: $NS (expected 135)"

echo
echo "[$(date)] === Step 2: annotate_vcf (THIS IS THE BIG ONE — ~30 min, ~20 GB RAM) ==="
if [ ! -s ${OUTPREFIX}.vcf.done ]; then
  $PY -u $ANNOTATE -vcf $VCF_CHR1 -gfa $GFA -o $OUTPREFIX 2>&1 | tail -10
  touch ${OUTPREFIX}.vcf.done
fi
echo "  annotated multi-allelic: $(grep -vc '^#' ${OUTPREFIX}.vcf) records"
echo "  biallelic catalog:       $(grep -vc '^#' ${OUTPREFIX}_biallelic.vcf) records"

echo
echo "[$(date)] === Step 3: sort + bgzip + index ==="
for f in ${OUTPREFIX} ${OUTPREFIX}_biallelic; do
  if [ ! -s ${f}.sorted.vcf.gz ]; then
    awk '$1 ~ /^#/ {print $0; next} {print $0 | "sort -k1,1 -k2,2n"}' ${f}.vcf > ${f}.sorted.vcf
    $BGZIP -f -c ${f}.sorted.vcf > ${f}.sorted.vcf.gz
    $TABIX -p vcf ${f}.sorted.vcf.gz
  fi
done
ls -lh ${OUTPREFIX}.sorted.vcf.gz ${OUTPREFIX}_biallelic.sorted.vcf.gz

echo
echo "[$(date)] DONE A1 (annotate Chr1)"
