#!/bin/bash
#SBATCH --job-name=archB_merge
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/jobB_merge_%j.out
#SBATCH --error=logs/jobB_merge_%j.err
mkdir -p logs
set -euo pipefail

# Job B v2: convert cactus_78 to biallelic using 135-catalog, merge with PG_153 from Job A,
# post-merge AN=0 filter, panel-level AC sanity.
# FIX vs v1: don't try to subset 135-annotated VCF by 1001G IDs (sample names don't match
# because 135-VCF has original cactus assembly names). Instead, use cactus_78.vcf.gz
# (already in 1001G IDs) and transfer INFO/ID from the 135-catalog.

cd /global/scratch/users/tbellg/kmate/scratch/arch3_test

BCF=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bcftools
BGZIP=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bgzip
TABIX=/global/home/users/tbellg/miniforge3/envs/gwas/bin/tabix
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CONVERT=/global/scratch/users/tbellg/kmate/external_tools/pangenie-tools/pipelines/run-from-callset/scripts/convert-to-biallelic.py
TRANSFER=/global/scratch/users/tbellg/kmate/scratch/arch3_test/transfer_id_annotation.py

CACTUS78_RAW=/global/scratch/users/tbellg/kmate/panel/pangenie_genotyping/data/v3qc_tmp/cactus_78.vcf.gz
CACTUS_ANNOT_135=full135_test_annotated.sorted.vcf.gz
BIAL_CATALOG=full135_test_annotated_biallelic.sorted.vcf.gz
PG_HAP_FROM_A=pg_153_test_haploid.vcf.gz

[ -s $PG_HAP_FROM_A ] || { echo "ERROR: missing PG haploid from Job A"; exit 1; }
[ -s $CACTUS_ANNOT_135 ] || { echo "ERROR: missing 135-annotated cactus"; exit 1; }

echo "[$(date)] === Step 1: subset cactus_78 to test region ==="
CACTUS78_TEST=cactus_78_test_region.vcf.gz
if [ ! -s $CACTUS78_TEST ]; then
  $BCF view -r Chr1:5800000-14000000 $CACTUS78_RAW -Oz -o $CACTUS78_TEST
  $TABIX -p vcf $CACTUS78_TEST
fi
N=$($BCF view -H $CACTUS78_TEST | wc -l)
NS=$($BCF query -l $CACTUS78_TEST | wc -l)
echo "  cactus_78 test region records: $N, samples: $NS"

echo
echo "[$(date)] === Step 2: transfer INFO/ID from 135-catalog ==="
CACTUS78_ANNOT=cactus_78_test_annotated.vcf
$PY $TRANSFER --cactus $CACTUS_ANNOT_135 --pg $CACTUS78_TEST --out $CACTUS78_ANNOT

echo
echo "[$(date)] === Step 3: convert-to-biallelic on cactus_78 using 135-catalog ==="
CACTUS78_BIAL=cactus_78_test_biallelic.vcf
cat $CACTUS78_ANNOT | $PY $CONVERT $BIAL_CATALOG > $CACTUS78_BIAL 2> cactus_convert_v2.log
echo "  exit: $?"
echo "  biallelic records: $(grep -vc '^#' $CACTUS78_BIAL)"

echo
echo "[$(date)] === Step 4: sort + bgzip + index cactus_78 biallelic ==="
CACTUS78_HAP=cactus_78_test_haploid.vcf.gz
awk '$1 ~ /^#/ {print $0; next} {print $0 | "sort -k1,1 -k2,2n"}' $CACTUS78_BIAL > cactus_78_test_biallelic.sorted.vcf
$BGZIP -f -c cactus_78_test_biallelic.sorted.vcf > $CACTUS78_HAP
$TABIX -p vcf $CACTUS78_HAP
echo "  $CACTUS78_HAP: $(ls -la $CACTUS78_HAP | awk '{print $5}') bytes"

echo
echo "[$(date)] === Step 5: bcftools merge (both haploid biallelic) ==="
MERGED=merged_88_test_v2.vcf.gz
$BCF merge $CACTUS78_HAP $PG_HAP_FROM_A --threads 4 -Oz -o $MERGED
$TABIX -p vcf $MERGED
N_MERGED=$($BCF view -H $MERGED | wc -l)
N_SAMP=$($BCF query -l $MERGED | wc -l)
echo "  merged records: $N_MERGED"
echo "  merged samples: $N_SAMP (expected 78 + 153 = 231)"

echo
echo "[$(date)] === Step 6: fill-tags on merged ==="
MERGED_FILLED=merged_231_test_filled.vcf.gz
$BCF +fill-tags $MERGED --threads 4 -Oz -o $MERGED_FILLED -- -t AC,AN,F_MISSING
$TABIX -p vcf $MERGED_FILLED

echo
echo "[$(date)] === Step 7: post-merge AN=0 filter ==="
MERGED_FINAL=merged_231_test_final.vcf.gz
N_PRE=$($BCF view -H $MERGED_FILLED | wc -l)
$BCF view -e 'INFO/AN=0' $MERGED_FILLED --threads 4 -Oz -o $MERGED_FINAL
$TABIX -p vcf $MERGED_FINAL
N_POST=$($BCF view -H $MERGED_FINAL | wc -l)
echo "  records pre-AN=0:  $N_PRE"
echo "  records post-AN=0: $N_POST"
echo "  dropped: $((N_PRE - N_POST))"

echo
echo "[$(date)] === Step 8: panel-level AC at 3 spot-checks ==="
echo "Expected (xwu reference):"
echo "  5870018 T->A: AC ~5/231 (~2%)"
echo "  10421645 T->C: AC ~229/231 (~99%)"
echo
$PY <<'PYEOF'
import gzip
SPOT = [(5870018, 'T', 'A'), (10421645, 'T', 'C'),
        (13843898, 'C', 'T'), (13843898, 'CT', 'TC'), (13843898, 'CT', 'TG')]
results = {}
samples = None
with gzip.open('merged_231_test_final.vcf.gz', 'rt') as f:
    for line in f:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            samples = line.rstrip().split('\t')[9:]
            continue
        parts = line.rstrip().split('\t')
        key = (int(parts[1]), parts[3], parts[4])
        if key in SPOT:
            results[key] = parts
print(f"merged panel samples: {len(samples)}")
for key in SPOT:
    if key not in results:
        print(f"  Chr1:{key[0]} {key[1]}>{key[2]}: not present")
        continue
    parts = results[key]
    gts = parts[9:]
    ac = sum(1 for g in gts if g == '1')
    an = sum(1 for g in gts if g in ('0','1'))
    af = ac/an if an>0 else 0
    print(f"  Chr1:{key[0]} {key[1]}>{key[2]}: panel_AC={ac}/{an}  AF={af:.4f}")

print()
print("=== OR-merged at coord 13843898 (sum carriers across 3 records with T at coord) ===")
target_pos = 13843898
carriers = [False]*len(samples)
called = [False]*len(samples)
records_seen = []
with gzip.open('merged_231_test_final.vcf.gz', 'rt') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.rstrip().split('\t')
        if int(parts[1]) != target_pos: continue
        if not parts[4] or parts[4][0] != 'T': continue
        records_seen.append((parts[3], parts[4]))
        for i, gt in enumerate(parts[9:]):
            if gt == '1': carriers[i] = True; called[i] = True
            elif gt == '0': called[i] = True
print(f"  records contributing: {records_seen}")
print(f"  OR-merged carriers: {sum(carriers)}/{sum(called)} called  (total {len(samples)})")
PYEOF

echo
echo "[$(date)] DONE Job B v2"
ls -lh merged_231_test_final.vcf.gz cactus_78_test_haploid.vcf.gz
