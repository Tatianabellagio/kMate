#!/bin/bash
#SBATCH --job-name=pg_miss2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/pg_miss2_%j.out
#SBATCH --error=logs/pg_miss2_%j.err

mkdir -p logs
set -uo pipefail
BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping/data/v3qc_v2
OUT=$BASE/missingness_report

PG=$BASE/pangenie_153_qc_v2.vcf.gz
CACTUS=$BASE/cactus_78_bi.vcf.gz
MERGED=$BASE/founders_231_v3qc_v2.vcf.gz

# ---- Per-sample missing rate (PG) — write stats once, parse twice ----
echo "[$(date)] Section 3: per-sample missing rate on PG"
$BCF stats -s- $PG > $OUT/pg_bcfstats.txt
echo "  bcftools stats done"

echo ""
echo "Header: sample  nRefHom  nNonRefHom  nHets  nMissing"
awk '/^PSC/{
  sample=$3; nrefhom=$4; nalthom=$5; nhet=$6; nmissing=$13
  total = nrefhom + nalthom + nhet + nmissing
  if (total > 0) printf "%s\t%d\t%d\t%.4f\n", sample, nmissing, total, nmissing/total
}' $OUT/pg_bcfstats.txt > $OUT/pg_per_sample_missing.tsv

awk -F'\t' '{
  v = $4 + 0
  sum += v; n++
  if (n == 1 || v < min) min = v
  if (v > max) max = v
} END {
  printf "Samples: %d  mean_miss: %.4f  min: %.4f  max: %.4f\n", n, sum/n, min, max
}' $OUT/pg_per_sample_missing.tsv

echo ""
echo "Top 10 highest-missing PG samples:"
sort -t$'\t' -k4 -g -r $OUT/pg_per_sample_missing.tsv | head -10 | awk -F'\t' '{printf "  %-20s miss=%.4f  (%d / %d)\n", $1, $4, $2, $3}'

echo ""
echo "Bottom 10 least-missing PG samples:"
sort -t$'\t' -k4 -g $OUT/pg_per_sample_missing.tsv | head -10 | awk -F'\t' '{printf "  %-20s miss=%.4f  (%d / %d)\n", $1, $4, $2, $3}'

# ---- Per-chrom missing rate ----
echo ""
echo "[$(date)] Section 4: per-chrom F_MISSING from pg_class.tsv"
awk -F'\t' '$5 != "." {
  v = $5 + 0
  sum[$1] += v; n[$1]++
} END {
  printf "  %-8s %12s %12s\n", "chrom", "n_records", "mean_FMISS"
  for (c in n) printf "  %-8s %12d %12.4f\n", c, n[c], sum[c]/n[c]
}' $OUT/pg_class.tsv | sort

# ---- Cactus control ----
echo ""
echo "[$(date)] Section 5: cactus_78_bi F_MISSING (control — should be ~0)"
$BCF query -f '%INFO/F_MISSING\n' $CACTUS > $OUT/cactus_fmissing.txt
echo "  Total: $(wc -l < $OUT/cactus_fmissing.txt) records"
awk '$1 != "." {v=$1+0; sum+=v; n++; if (v>max) max=v; if (v==0) zero++} END {printf "  mean: %.6f, max: %.4f, %%zero: %.2f%%\n", sum/n, max, 100*zero/n}' $OUT/cactus_fmissing.txt

# ---- Merged 231 panel ----
echo ""
echo "[$(date)] Section 6: founders_231_v3qc_v2 (231 samples, what cn_full sees)"
$BCF query -f '%INFO/F_MISSING\n' $MERGED > $OUT/merged_fmissing.txt
echo "  Total: $(wc -l < $OUT/merged_fmissing.txt) records"
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
  printf "    0.00     : %d (%.2f%%)\n", bin0, 100*bin0/n
  printf "    0.00-0.01: %d (%.2f%%)\n", bin1, 100*bin1/n
  printf "    0.01-0.05: %d (%.2f%%)\n", bin2, 100*bin2/n
  printf "    0.05-0.10: %d (%.2f%%)\n", bin3, 100*bin3/n
  printf "    0.10-0.20: %d (%.2f%%)\n", bin4, 100*bin4/n
  printf "    0.20-0.50: %d (%.2f%%)\n", bin5, 100*bin5/n
  printf "    0.50-1.00: %d (%.2f%%)\n", bin6, 100*bin6/n
}' $OUT/merged_fmissing.txt

# ---- Bonus: total missing CELL count broken by class ----
echo ""
echo "[$(date)] Section 7: total missing-cell budget by class (153 samples per record)"
awk -F'\t' '$5 != "." {
  ref_len = length($3); alt_len = length($4)
  if (ref_len == 1 && alt_len == 1) class = "SNP"
  else if (ref_len == alt_len) class = "MNP"
  else {
    diff = (ref_len > alt_len) ? ref_len - alt_len : alt_len - ref_len
    class = (diff >= 50) ? "SV>=50bp" : "indel<50bp"
  }
  v = $5 + 0
  count[class]++
  cells_total[class] += 153
  cells_missing[class] += v * 153
} END {
  printf "  %-12s %12s %15s %18s %10s\n", "class", "n_records", "total_cells", "missing_cells", "%miss"
  for (c in count) {
    printf "  %-12s %12d %15d %18d %10.3f\n", c, count[c], cells_total[c], int(cells_missing[c]), 100*cells_missing[c]/cells_total[c]
  }
}' $OUT/pg_class.tsv

echo ""
echo "[$(date)] DONE"
