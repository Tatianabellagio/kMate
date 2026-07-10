#!/bin/bash
# Fast merged-VCF stats via bcftools query streaming.
# Per-record (carriers, AC, AN, F_MISSING) by full panel + each cohort.
# Per-sample: streamed counter aggregation.
set -eo pipefail

BASE=/global/scratch/users/tbellg/hapfire_sv
VCF=${VCF:-$BASE/pangenie_genotyping/data/merged/founders_231_chr.vcf.gz}
OUT=${OUT:-$BASE/preprocess_qc/output/merged_stats}
mkdir -p $OUT

BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

# Build cohort sample lists
ALL_SAMPLES=$($BCF query -l $VCF)
CACTUS_OVERLAP=$BASE/imputation/work/sample_rename.txt
awk '{print $2}' $CACTUS_OVERLAP | sort -u > $OUT/cactus_80.txt
echo "$ALL_SAMPLES" | sort -u > $OUT/all_231.txt
comm -23 $OUT/all_231.txt $OUT/cactus_80.txt > $OUT/pangenie_151.txt
echo "Sample lists: cactus_80=$(wc -l < $OUT/cactus_80.txt)  pangenie_151=$(wc -l < $OUT/pangenie_151.txt)  all_231=$(wc -l < $OUT/all_231.txt)"

# Per-record TSV via bcftools query — outputs (chrom, pos, ref_len, alt_max_len, GT-string).
# Then a single Python pass aggregates: per-record AC/AN/missing for full panel + each cohort,
# plus per-sample tallies. Single 5.21M-record stream.

echo "[$(date)] streaming bcftools query -> python aggregator"
time $BCF query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' $VCF | \
    OUT_PREFIX=${OUT_PREFIX:-$OUT/founders_231}
    $PY $BASE/preprocess_qc/scripts/merged_vcf_stats_aggregate.py \
        --cactus-samples $OUT/cactus_80.txt \
        --pangenie-samples $OUT/pangenie_151.txt \
        --out-prefix $OUT_PREFIX \
        --vcf $VCF

echo "[$(date)] DONE"
ls -lh $OUT/founders_231*
