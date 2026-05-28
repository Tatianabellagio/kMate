#!/bin/bash
#SBATCH --job-name=arch3_pg_v3
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/pg_v3_%j.out
#SBATCH --error=logs/pg_v3_%j.err
mkdir -p logs
set -euo pipefail

cd /global/scratch/users/tbellg/hapfire_sv/scratch/arch3_test

BCF=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bcftools
BGZIP=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/gwas/bin/tabix
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CONVERT=/global/scratch/users/tbellg/hapfire_sv/external_tools/pangenie-tools/pipelines/run-from-callset/scripts/convert-to-biallelic.py
TRANSFER=/global/scratch/users/tbellg/hapfire_sv/scratch/arch3_test/transfer_id_annotation.py

CACTUS_ANNOT=full135_test_annotated.sorted.vcf.gz   # from job 63045
BIAL_CATALOG=full135_test_annotated_biallelic.sorted.vcf.gz
PG_TEST=pg_100001_test_chr1_5_14M.vcf.gz

[ -s $CACTUS_ANNOT ] || { echo "ERROR: missing $CACTUS_ANNOT (run job 63045 first)"; exit 1; }
[ -s $PG_TEST ] || { echo "ERROR: missing $PG_TEST"; exit 1; }

echo "[$(date)] === Step 1: custom Python ID transfer ==="
PG_ANNOT=pg_100001_test_annotated_v3.vcf
$PY $TRANSFER --cactus $CACTUS_ANNOT --pg $PG_TEST --out $PG_ANNOT

echo
echo "[$(date)] === Step 2: spot-check INFO/ID is now populated ==="
echo "Records with non-empty INFO/ID:"
awk -F'\t' '
  !/^#/ {
    n_total++;
    if ($8 ~ /(^|;)ID=[^;]/) n_with_id++;
  }
  END { print "  "n_with_id" / "n_total" ("100*n_with_id/n_total"%)" }
' $PG_ANNOT

echo
echo "Spot-check Chr1:5870018 (over-counting position):"
awk -F'\t' '!/^#/ && $2==5869846' $PG_ANNOT | head -1 | awk -F'\t' '{
  print "  REF_len=" length($4)
  for (n=split($8,a,";"); i++<n;) if (a[i] ~ /^ID=/) print "  INFO/ID first 250 chars: " substr(a[i],4,250)
}'

echo
echo "[$(date)] === Step 3: convert-to-biallelic ==="
PG_BIAL=pg_100001_per_sample_biallelic_v3.vcf
cat $PG_ANNOT | $PY $CONVERT $BIAL_CATALOG > $PG_BIAL 2> pg_convert_v3.log
echo "convert exit: $?"
echo "PG biallelic records: $(grep -vc '^#' $PG_BIAL)"
tail -5 pg_convert_v3.log

echo
echo "[$(date)] === Step 4: spot-check GTs at 3 positions ==="
for entry in "5870018:T:A" "10421645:T:C" "13843898:C:T" "13843898:CT:TC" "13843898:CT:TG"; do
  IFS=':' read -r pos ref alt <<< "$entry"
  echo "--- Chr1:$pos REF=$ref ALT=$alt ---"
  awk -v p=$pos -v r=$ref -v a=$alt '
    !/^#/ && $2==p && $4==r && $5==a {
      gt=$10;
      split(gt, fmt, ":"); g=fmt[1];
      print "  sample 100001 GT="g
    }
  ' $PG_BIAL
done

echo
echo "[$(date)] DONE"
