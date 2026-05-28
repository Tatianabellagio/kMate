#!/bin/bash
#SBATCH --job-name=archD_centro
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=logs/jobD_centro_%j.out
#SBATCH --error=logs/jobD_centro_%j.err
mkdir -p logs
set -euo pipefail

# Job D v2: full Arch 3 pipeline on centromere region Chr1:14M-17M.
# FIX vs v1: for cactus side, use cactus_78.vcf.gz (1001G IDs) + transfer_id from centromere catalog
# instead of bcftools view -S on the 135-VCF (which has different sample names).

cd /global/scratch/users/tbellg/kmate/scratch/arch3_test

BCF=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bcftools
BGZIP=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/gwas/bin/tabix
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
ANNOTATE=/global/scratch/users/tbellg/kmate/external/genotyping-pipelines/prepare-vcf-MC/workflow/scripts/annotate_vcf.py
CONVERT=/global/scratch/users/tbellg/kmate/external/pangenie-tools/pipelines/run-from-callset/scripts/convert-to-biallelic.py
TRANSFER=/global/scratch/users/tbellg/kmate/scratch/arch3_test/transfer_id_annotation.py

GFA=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.gfa.gz
FULL_135_VCF=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.vcf.gz
PG_RAW=/global/scratch/users/tbellg/kmate/panel/pangenie_genotyping/data/v3qc_tmp/pangenie_153_raw.vcf.gz
CACTUS78_RAW=/global/scratch/users/tbellg/kmate/panel/pangenie_genotyping/data/v3qc_tmp/cactus_78.vcf.gz

REGION="Chr1:14000000-17000000"
OUTDIR=centromere
mkdir -p $OUTDIR
cd $OUTDIR

# Reuse Step 1-3 outputs from previous run (centro_135_annotated.sorted.vcf.gz exists)
if [ -s centro_135_annotated.sorted.vcf.gz ] && [ -s centro_135_annotated_biallelic.sorted.vcf.gz ]; then
  echo "[$(date)] === Reusing annotated catalog from prior run ==="
  ls -lh centro_135_annotated.sorted.vcf.gz centro_135_annotated_biallelic.sorted.vcf.gz
else
  echo "[$(date)] === Step 1: subset 135-sample VCF to centromere ==="
  $BCF view -r $REGION $FULL_135_VCF > centro_135.vcf
  echo "  records: $($BCF view -H centro_135.vcf | wc -l)"

  echo
  echo "[$(date)] === Step 2: annotate_vcf on centromere ==="
  $PY -u $ANNOTATE -vcf centro_135.vcf -gfa $GFA -o centro_135_annotated 2>&1 | tail -5

  echo
  echo "[$(date)] === Step 3: sort + bgzip + index ==="
  for f in centro_135_annotated centro_135_annotated_biallelic; do
    awk '$1 ~ /^#/ {print $0; next} {print $0 | "sort -k1,1 -k2,2n"}' ${f}.vcf > ${f}.sorted.vcf
    $BGZIP -f -c ${f}.sorted.vcf > ${f}.sorted.vcf.gz
    $TABIX -p vcf ${f}.sorted.vcf.gz
  done
fi

# Reuse PG haploid from prior run if it exists
if [ -s pg_153_centro_haploid.vcf.gz ]; then
  echo "[$(date)] === Reusing PG haploid from prior run ==="
  ls -lh pg_153_centro_haploid.vcf.gz
else
  echo "[$(date)] === Step 4: PG side — subset, transfer ID, convert, fill-tags, V4 filter, haploidize ==="
  $BCF view -r $REGION $PG_RAW -Oz -o pg_153_centro.vcf.gz
  $TABIX -p vcf pg_153_centro.vcf.gz

  $PY $TRANSFER --cactus centro_135_annotated.sorted.vcf.gz --pg pg_153_centro.vcf.gz --out pg_153_centro_annotated.vcf

  $BCF view pg_153_centro_annotated.vcf 2>/dev/null | $PY $CONVERT centro_135_annotated_biallelic.sorted.vcf.gz > pg_153_centro_biallelic.vcf 2> pg_convert.log

  $BCF +fill-tags pg_153_centro_biallelic.vcf --threads 4 -Oz -o pg_153_centro_filled.vcf.gz -- -t AC,AN,AC_Het,F_MISSING

  # -------------------------------------------------------------------------
  # Haploidize PG side — het handling design decision (see jobA_pg_full_pipeline.sh
  # for full rationale). Summary:
  #   0/0 → 0 (REF call preserved)
  #   1/1 → 1 (ALT carrier preserved)
  #   ./. → . (missing preserved)
  #   0/1 → . (HET → MISSING; honest, follows Arouisse 2020 Arabidopsis precedent)
  #   else → . (defensive)
  # NO V4 row filter applied. All records preserved; per-cell het→missing is
  # the only PG-side het handling.
  # -------------------------------------------------------------------------
  $BCF view --threads 2 -Ov pg_153_centro_filled.vcf.gz | \
  awk 'BEGIN{OFS="\t"}
      /^##/ { print; next }
      /^#CHROM/ { print; next }
      {
          $9 = "GT"
          for (i=10; i<=NF; i++) {
              split($i, parts, ":"); g = parts[1]
              if (g == "0/0" || g == "0|0") { $i = "0" }
              else if (g == "1/1" || g == "1|1") { $i = "1" }
              else if (g == "./." || g == ".|." || g == ".") { $i = "." }
              else if (g ~ /^0[\/|]1$/ || g ~ /^1[\/|]0$/) { $i = "." }   # het → MISSING (Arouisse 2020)
              else { $i = "." }
          }
          print
      }' | $BGZIP -@ 4 -c > pg_153_centro_haploid.vcf.gz
  $TABIX -p vcf pg_153_centro_haploid.vcf.gz
fi

echo
echo "[$(date)] === Step 5: cactus side — subset cactus_78.vcf.gz directly (NOT bcftools view -S on 135-VCF) ==="
# Subset cactus_78.vcf.gz (which is already in 1001G IDs) to centromere
$BCF view -r $REGION $CACTUS78_RAW -Oz -o cactus_78_centro_region.vcf.gz
$TABIX -p vcf cactus_78_centro_region.vcf.gz
N=$($BCF view -H cactus_78_centro_region.vcf.gz | wc -l)
NS=$($BCF query -l cactus_78_centro_region.vcf.gz | wc -l)
echo "  cactus_78 centromere records: $N, samples: $NS"

# Transfer INFO/ID from centromere catalog to cactus_78
CACTUS78_CENTRO_ANNOT=cactus_78_centro_annotated.vcf
$PY $TRANSFER --cactus centro_135_annotated.sorted.vcf.gz --pg cactus_78_centro_region.vcf.gz --out $CACTUS78_CENTRO_ANNOT

# Convert to biallelic
cat $CACTUS78_CENTRO_ANNOT | $PY $CONVERT centro_135_annotated_biallelic.sorted.vcf.gz > centro_cactus78_biallelic.vcf 2> cactus_convert.log
echo "  cactus biallelic records: $(grep -vc '^#' centro_cactus78_biallelic.vcf)"

awk '$1 ~ /^#/ {print $0; next} {print $0 | "sort -k1,1 -k2,2n"}' centro_cactus78_biallelic.vcf > centro_cactus78_biallelic.sorted.vcf
$BGZIP -f -c centro_cactus78_biallelic.sorted.vcf > centro_cactus78_haploid.vcf.gz
$TABIX -p vcf centro_cactus78_haploid.vcf.gz
ls -lh centro_cactus78_haploid.vcf.gz

echo
echo "[$(date)] === Step 6: merge + filter + summary ==="
$BCF merge centro_cactus78_haploid.vcf.gz pg_153_centro_haploid.vcf.gz --threads 4 -Oz -o centro_merged_231.vcf.gz
$TABIX -p vcf centro_merged_231.vcf.gz
$BCF +fill-tags centro_merged_231.vcf.gz --threads 4 -Oz -o centro_merged_231_filled.vcf.gz -- -t AC,AN,F_MISSING
$TABIX -p vcf centro_merged_231_filled.vcf.gz
N_PRE=$($BCF view -H centro_merged_231_filled.vcf.gz | wc -l)
$BCF view -e 'INFO/AN=0' centro_merged_231_filled.vcf.gz --threads 4 -Oz -o centro_merged_231_final.vcf.gz
N_POST=$($BCF view -H centro_merged_231_final.vcf.gz | wc -l)
echo "  centromere merged records: $N_PRE -> $N_POST after AN=0 filter"
echo "  centromere samples: $($BCF query -l centro_merged_231_final.vcf.gz | wc -l) (expected 231)"

echo
echo "[$(date)] === Step 7: F_MISSING distribution on centromere final ==="
$BCF query -f '%INFO/F_MISSING\n' centro_merged_231_final.vcf.gz | $PY -c "
import sys, numpy as np
fm = np.array([float(x.strip()) for x in sys.stdin if x.strip() and x.strip() != '.'])
print(f'  N records: {len(fm):,}')
for lo, hi in [(0,0.01), (0.01,0.05), (0.05,0.1), (0.1,0.3), (0.3,1.01)]:
    n = ((fm >= lo) & (fm < hi)).sum()
    print(f'  F_MISSING in [{lo:.2f}, {hi:.2f}): {n:>7,} ({100*n/len(fm):.1f}%)')"

echo
echo "[$(date)] DONE Job D v2"
