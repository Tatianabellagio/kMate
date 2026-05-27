#!/bin/bash
#SBATCH --job-name=pg_missing
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/pg_missing_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/pg_missing_%j.err

# Quantify missingness in pangenie_153_qc_v2 across:
#  - global F_MISSING distribution
#  - per-variant-class (SNP / indel<50bp / SV>=50bp / MNP)
#  - per-sample missing rate
#  - per-chrom
# Compare against cactus_78_bi (should be ~0% missing).
#
# Purpose: decide if missingness is meaningful enough to justify Beagle imputation.
# If PG missing rate <5%, then k-mer imbalance comes mostly from long-read vs
# short-read modality (not from ./.), and Beagle won't fix the bias.

set -euo pipefail
BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/data/v3qc_v2
OUT=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/data/v3qc_v2/missingness_report

PG=$BASE/pangenie_153_qc_v2.vcf.gz
CACTUS=$BASE/cactus_78_bi.vcf.gz
MERGED=$BASE/founders_231_v3qc_v2.vcf.gz

mkdir -p $OUT

echo "=========================================================="
echo "[$(date)] PG missingness analysis"
echo "=========================================================="

# ---- 1. Global F_MISSING distribution on PG (already filled-in v3qc_v2) ----
echo ""
echo "--- 1. PG global F_MISSING distribution (153 samples) ---"
$BCF query -f '%INFO/F_MISSING\n' $PG > $OUT/pg_fmissing.txt
echo "Total records: $(wc -l < $OUT/pg_fmissing.txt)"
echo "Histogram of F_MISSING:"
awk '{
  if ($1 == ".") next
  v = $1 + 0
  if (v == 0) bin = "0.00"
  else if (v < 0.01) bin = "0.00-0.01"
  else if (v < 0.05) bin = "0.01-0.05"
  else if (v < 0.10) bin = "0.05-0.10"
  else if (v < 0.20) bin = "0.10-0.20"
  else if (v < 0.30) bin = "0.20-0.30"
  else if (v < 0.50) bin = "0.30-0.50"
  else bin = "0.50-1.00"
  count[bin]++
  total++
} END {
  for (b in count) printf "  %-12s %10d (%.2f%%)\n", b, count[b], 100*count[b]/total
}' $OUT/pg_fmissing.txt | sort

echo ""
echo "Summary stats of F_MISSING (PG):"
awk '$1 != "." {v=$1+0; sum+=v; n++; if (v>max) max=v} END {printf "  mean: %.4f, max: %.4f, n: %d\n", sum/n, max, n}' $OUT/pg_fmissing.txt

# ---- 2. By variant class ----
echo ""
echo "--- 2. PG F_MISSING by variant class ---"
$BCF query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/F_MISSING\n' $PG > $OUT/pg_class.tsv
awk -F'\t' '$5 != "." {
  ref_len = length($3)
  alt_len = length($4)
  if (ref_len == 1 && alt_len == 1) class = "SNP"
  else if (ref_len == alt_len) class = "MNP"
  else {
    diff = (ref_len > alt_len) ? ref_len - alt_len : alt_len - ref_len
    class = (diff >= 50) ? "SV>=50bp" : "indel<50bp"
  }
  v = $5 + 0
  count[class]++
  sum[class] += v
  if (v == 0) zero_miss[class]++
  if (v > 0) any_miss[class]++
} END {
  printf "  %-12s %12s %12s %12s %12s\n", "class", "n", "mean_FMISS", "%zero", "%any_miss"
  for (c in count) {
    printf "  %-12s %12d %12.4f %12.2f %12.2f\n", c, count[c], sum[c]/count[c], 100*zero_miss[c]/count[c], 100*any_miss[c]/count[c]
  }
}' $OUT/pg_class.tsv

# ---- 3. Per-sample missing rate (PG only) ----
echo ""
echo "--- 3. Per-sample missing rate (PG) ---"
$BCF stats -s- $PG | grep '^PSC' | head -1
$BCF stats -s- $PG | awk '/^PSC/{
  sample=$3; nrefhom=$4; nalthom=$5; nhet=$6; ntsv=$7; nindel=$8; avgdp=$9; nsingle=$10; nhapref=$11; nhapalt=$12; nmissing=$13
  total = nrefhom + nalthom + nhet + nmissing
  if (total > 0) printf "%s\t%d\t%d\t%.4f\n", sample, nmissing, total, nmissing/total
}' > $OUT/pg_per_sample_missing.tsv

echo "Per-sample missing rate distribution:"
awk -F'\t' '{
  v = $4 + 0
  sum += v; n++
  if (n == 1 || v < min) min = v
  if (v > max) max = v
  arr[NR] = v
} END {
  printf "  samples: %d, mean: %.4f, min: %.4f, max: %.4f\n", n, sum/n, min, max
}' $OUT/pg_per_sample_missing.tsv

echo ""
echo "Top 10 highest-missing samples:"
sort -t$'\t' -k4 -g -r $OUT/pg_per_sample_missing.tsv | head -10 | awk -F'\t' '{printf "  %-20s miss=%.4f  (%d / %d)\n", $1, $4, $2, $3}'

echo ""
echo "Bottom 10 (least-missing samples):"
sort -t$'\t' -k4 -g $OUT/pg_per_sample_missing.tsv | head -10 | awk -F'\t' '{printf "  %-20s miss=%.4f  (%d / %d)\n", $1, $4, $2, $3}'

# ---- 4. Per-chrom missing rate ----
echo ""
echo "--- 4. PG F_MISSING per chrom ---"
awk -F'\t' '$5 != "." {
  v = $5 + 0
  sum[$1] += v; n[$1]++
} END {
  printf "  %-8s %12s %12s\n", "chrom", "n_records", "mean_FMISS"
  for (c in n) printf "  %-8s %12d %12.4f\n", c, n[c], sum[c]/n[c]
}' $OUT/pg_class.tsv | sort

# ---- 5. Cactus comparison (control) ----
echo ""
echo "--- 5. Cactus F_MISSING distribution (control — should be ~0) ---"
$BCF query -f '%INFO/F_MISSING\n' $CACTUS > $OUT/cactus_fmissing.txt
echo "Total records: $(wc -l < $OUT/cactus_fmissing.txt)"
awk '$1 != "." {v=$1+0; sum+=v; n++; if (v>max) max=v; if (v==0) zero++} END {printf "  mean: %.6f, max: %.4f, %%zero: %.2f%%\n", sum/n, max, 100*zero/n}' $OUT/cactus_fmissing.txt

# ---- 6. Cross-check on the 231-founder merged VCF ----
echo ""
echo "--- 6. founders_231_v3qc_v2 (231 samples) F_MISSING — what cn_full actually sees ---"
$BCF query -f '%INFO/F_MISSING\n' $MERGED > $OUT/merged_fmissing.txt
awk '$1 != "." {
  v=$1+0; sum+=v; n++
  if (v==0) bin0++
  else if (v < 0.01) bin1++
  else if (v < 0.05) bin2++
  else if (v < 0.10) bin3++
  else if (v < 0.20) bin4++
  else if (v < 0.50) bin5++
  else bin6++
} END {
  printf "  total: %d, mean F_MISSING: %.4f\n", n, sum/n
  printf "    0.00    : %d (%.2f%%)\n", bin0, 100*bin0/n
  printf "    0.00-0.01: %d (%.2f%%)\n", bin1, 100*bin1/n
  printf "    0.01-0.05: %d (%.2f%%)\n", bin2, 100*bin2/n
  printf "    0.05-0.10: %d (%.2f%%)\n", bin3, 100*bin3/n
  printf "    0.10-0.20: %d (%.2f%%)\n", bin4, 100*bin4/n
  printf "    0.20-0.50: %d (%.2f%%)\n", bin5, 100*bin5/n
  printf "    0.50-1.00: %d (%.2f%%)\n", bin6, 100*bin6/n
}' $OUT/merged_fmissing.txt

echo ""
echo "[$(date)] DONE"
ls -lh $OUT/
