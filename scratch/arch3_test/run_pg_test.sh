#!/bin/bash
#SBATCH --job-name=arch3_pg
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/pg_test_%j.out
#SBATCH --error=logs/pg_test_%j.err
mkdir -p logs
set -euo pipefail

cd /global/scratch/users/tbellg/kmate/scratch/arch3_test

BCF=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bcftools
BGZIP=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/gwas/bin/tabix
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CONVERT=/global/scratch/users/tbellg/kmate/external/pangenie-tools/pipelines/run-from-callset/scripts/convert-to-biallelic.py

PG_RAW=/global/scratch/users/tbellg/kmate/panel/pangenie_genotyping/data/genotyped/100001_genotyping.vcf.gz
CACTUS_ANNOT=cactus_78_test_annotated.sorted.vcf
BIAL_CATALOG=cactus_78_test_annotated_biallelic.sorted.vcf.gz

echo "[$(date)] === Step 1: subset PG to test region ==="
PG_TEST=pg_100001_test_chr1_5_14M.vcf.gz
if [ ! -s $PG_TEST ]; then
  $BCF view -r Chr1:5800000-14000000 $PG_RAW -Oz -o $PG_TEST
  $TABIX -p vcf $PG_TEST
fi
$BCF view -H $PG_TEST | wc -l
echo "PG records in test region: ^"

echo
echo "[$(date)] === Step 2: bgzip+index the annotated cactus VCF (annotation source) ==="
CACTUS_ANNOT_GZ=cactus_78_test_annotated.sorted.vcf.gz
if [ ! -s $CACTUS_ANNOT_GZ ]; then
  $BGZIP -c $CACTUS_ANNOT > $CACTUS_ANNOT_GZ
  $TABIX -p vcf $CACTUS_ANNOT_GZ
fi
ls -la $CACTUS_ANNOT_GZ

echo
echo "[$(date)] === Step 3: transfer INFO/ID from cactus to PG ==="
PG_ANNOT=pg_100001_test_annotated.vcf.gz
$BCF annotate -a $CACTUS_ANNOT_GZ -c INFO/ID $PG_TEST -Oz -o $PG_ANNOT
$TABIX -p vcf $PG_ANNOT

echo
echo "=== Verify INFO/ID populated in PG-annotated ==="
echo "Sample records with non-empty INFO/ID:"
$BCF view -H $PG_ANNOT 2>/dev/null | awk -F'\t' '
  {
    info=$8;
    if (info ~ /(^|;)ID=[^;]/) n_with_id++;
    n_total++;
  }
  END { print "  with INFO/ID: "n_with_id" / "n_total }
'

echo
echo "[$(date)] === Step 4: run convert-to-biallelic on annotated PG ==="
PG_BIAL=pg_100001_per_sample_biallelic.vcf
$BCF view $PG_ANNOT 2>/dev/null | $PY $CONVERT $BIAL_CATALOG > $PG_BIAL 2> pg_convert.log
echo "Exit code: $?"
echo "Records emitted: $(grep -vc '^#' $PG_BIAL)"

echo
echo "[$(date)] === Step 5: spot-check PG biallelic GTs at 3 positions ==="
for entry in "5870018:T:A" "10421645:T:C" "13843898:C:T" "13843898:CT:TC" "13843898:CT:TG"; do
  IFS=':' read -r pos ref alt <<< "$entry"
  echo "--- Chr1:$pos REF=$ref ALT=$alt ---"
  awk -v p=$pos -v r=$ref -v a=$alt '
    !/^#/ && $2==p && $4==r && $5==a {
      n_rec++;
      for (i=10; i<=NF; i++) {
        split($i, fmt, ":"); gt=fmt[1];
        gsub(/[|\/]/, " ", gt);
        n_a=split(gt, alleles, " ");
        any_alt=0; any_called=0;
        for (j=1; j<=n_a; j++) {
          if (alleles[j] != ".") any_called=1;
          if (alleles[j] != "." && alleles[j] != "0") any_alt=1;
        }
        if (any_alt) n_carriers++;
        if (any_called) n_called++;
      }
      print "  records="n_rec"  carriers="n_carriers"  called="n_called"  (single-sample VCF — should be 0 or 1)";
      n_carriers=0; n_called=0;
    }
  ' $PG_BIAL
done

echo
echo "[$(date)] DONE"
