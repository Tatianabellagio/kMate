#!/bin/bash
#SBATCH --job-name=arch3_annot135
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/annot135_%j.out
#SBATCH --error=logs/annot135_%j.err
mkdir -p logs
set -euo pipefail

cd /global/scratch/users/tbellg/kmate/scratch/arch3_test

BCF=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bcftools
BGZIP=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/gwas/bin/tabix
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
ANNOTATE_SCRIPT=/global/scratch/users/tbellg/kmate/external_tools/genotyping-pipelines/prepare-vcf-MC/workflow/scripts/annotate_vcf.py

GFA=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.gfa.gz
FULL_135_VCF=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.vcf.gz

echo "[$(date)] === Step 1: subset 135-sample VCF to test region ==="
VCF_TEST=full135_test_chr1_5_14M.vcf
if [ ! -s ${VCF_TEST}.done ]; then
  $BCF view -r Chr1:5800000-14000000 $FULL_135_VCF > $VCF_TEST
  touch ${VCF_TEST}.done
fi
N_REC=$($BCF view -H $VCF_TEST | wc -l)
echo "Records in 135-sample test VCF: $N_REC"
echo "Samples: $($BCF view -h $VCF_TEST | grep '^#CHROM' | awk '{print NF-9}') (expected 135)"

echo
echo "[$(date)] === Step 2: run annotate_vcf on 135-sample VCF ==="
OUTPREFIX=full135_test_annotated
if [ ! -s ${OUTPREFIX}.vcf.done ]; then
  $PY -u $ANNOTATE_SCRIPT -vcf $VCF_TEST -gfa $GFA -o $OUTPREFIX 2>&1 | tail -5
  touch ${OUTPREFIX}.vcf.done
fi

echo
echo "Annotated multi-allelic: $(grep -vc '^#' ${OUTPREFIX}.vcf) records"
echo "Biallelic catalog: $(grep -vc '^#' ${OUTPREFIX}_biallelic.vcf) records"

echo
echo "[$(date)] === Step 3: sort + bgzip + index ==="
for f in ${OUTPREFIX} ${OUTPREFIX}_biallelic; do
  if [ ! -s ${f}.sorted.vcf.gz ]; then
    awk '$1 ~ /^#/ {print $0; next} {print $0 | "sort -k1,1 -k2,2n"}' ${f}.vcf > ${f}.sorted.vcf
    $BGZIP -f -c ${f}.sorted.vcf > ${f}.sorted.vcf.gz
    $TABIX -p vcf ${f}.sorted.vcf.gz
  fi
done
ls -la ${OUTPREFIX}.sorted.vcf.gz ${OUTPREFIX}_biallelic.sorted.vcf.gz

echo
echo "[$(date)] === Step 4: test bcftools annotate on PG (the failure case) ==="
PG_RAW=/global/scratch/users/tbellg/kmate/pangenie_genotyping/data/genotyped/100001_genotyping.vcf.gz
PG_TEST=pg_100001_test_chr1_5_14M.vcf.gz
if [ ! -s $PG_TEST ]; then
  $BCF view -r Chr1:5800000-14000000 $PG_RAW -Oz -o $PG_TEST
  $TABIX -p vcf $PG_TEST
fi

PG_ANNOT=pg_100001_test_annotated_v2.vcf.gz
$BCF annotate -a ${OUTPREFIX}.sorted.vcf.gz -c INFO/ID $PG_TEST -Oz -o $PG_ANNOT 2>&1 | head -10
ANNOT_STATUS=$?
echo "bcftools annotate exit: $ANNOT_STATUS"
$TABIX -p vcf $PG_ANNOT

echo
echo "Records with non-empty INFO/ID in PG-annotated:"
$BCF view -H $PG_ANNOT 2>/dev/null | awk -F'\t' '
  BEGIN { n_with_id=0; n_total=0 }
  { n_total++; if ($8 ~ /(^|;)ID=[^;]/) n_with_id++ }
  END { print "  "n_with_id" / "n_total" ("100*n_with_id/n_total"%)" }
'

echo
echo "[$(date)] === Step 5: convert-to-biallelic on annotated PG ==="
CONVERT=/global/scratch/users/tbellg/kmate/external_tools/pangenie-tools/pipelines/run-from-callset/scripts/convert-to-biallelic.py
PG_BIAL=pg_100001_per_sample_biallelic_v2.vcf
$BCF view $PG_ANNOT 2>/dev/null | $PY $CONVERT ${OUTPREFIX}_biallelic.sorted.vcf.gz > $PG_BIAL 2> pg_convert_v2.log
echo "convert exit: $?"
echo "PG biallelic records: $(grep -vc '^#' $PG_BIAL)"
echo "PG convert log tail:"
tail -5 pg_convert_v2.log

echo
echo "[$(date)] === Step 6: spot-check sample 100001 at 3 positions ==="
for entry in "5870018:T:A" "10421645:T:C" "13843898:C:T" "13843898:CT:TC" "13843898:CT:TG"; do
  IFS=':' read -r pos ref alt <<< "$entry"
  echo "--- Chr1:$pos REF=$ref ALT=$alt ---"
  awk -v p=$pos -v r=$ref -v a=$alt '
    !/^#/ && $2==p && $4==r && $5==a {
      gt=$10;  # single-sample PG: only column 10
      split(gt, fmt, ":"); g=fmt[1];
      print "  sample 100001 GT="g
    }
  ' $PG_BIAL
done

echo
echo "[$(date)] DONE"
