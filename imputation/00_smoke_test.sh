#!/bin/bash
# Small interactive smoke test: subset everything to Chr1:1-2Mb, run Beagle,
# verify imputation works end-to-end before committing to the full run.
#
# Setup: 10 random cactus founders → 9 ref + 1 target with SVs masked.
# If concordance > 0.7 on this tiny region with sparse data, the workflow is sound.

set -euo pipefail
mkdir -p /carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/smoke
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/smoke

BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
BGZIP=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bgzip
TABIX=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/tabix
BEAGLE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/beagle.jar

# Inputs
CACTUS_BIALLELIC=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_test/raw_vcfbub_lv0_diploid.biallelic.norm.vcf.gz
GRENE=/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_updatedVCF/greneNet_final_v1.1.recode.vcf

REGION_CACTUS="Chr1:1-2000000"
REGION_GRENE="1:1-2000000"

# Build sample-rename map (we already have this from main framework, but redo here for self-contained test)
PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python
$PYTHON <<'PYEOF' > sample_rename.txt
import pandas as pd
panel = pd.read_csv('/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/sv_panel_to_accession_id.tsv', sep='\t')
asm_to_1001g = dict(zip(panel.Assembly_ID.astype(str), panel.Accession_ID.astype(str)))
grene = set(open('/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/vcf_samples_231.txt').read().split())
import subprocess
res = subprocess.run('/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools query -l '
                     '/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_test/raw_vcfbub_lv0_diploid.biallelic.norm.vcf.gz',
                     shell=True, capture_output=True, text=True)
for asm in res.stdout.strip().split():
    g = asm_to_1001g.get(asm)
    if g and g in grene:
        print(f'{asm}\t{g}')
PYEOF
echo -e "Chr1\t1\nChr2\t2\nChr3\t3\nChr4\t4\nChr5\t5" > chrom_rename.txt

# Pick 10 cactus founders
shuf --random-source=<(yes 42) sample_rename.txt | head -10 > pick10.txt
TEST_SAMPLE=$(head -1 pick10.txt | awk '{print $2}')
REF_SAMPLES=$(tail -9 pick10.txt | awk '{print $2}')
echo "Test (leave-out) sample: $TEST_SAMPLE"
echo "Reference samples (9): $REF_SAMPLES"

echo ""
echo "[$(date)] STEP 1: subset cactus VCF to Chr1:1-2Mb, rename, filter to 10 samples + SVs"
$BCF view -r $REGION_CACTUS $CACTUS_BIALLELIC --threads 4 \
    | $BCF reheader --samples sample_rename.txt \
    | $BCF annotate --rename-chrs chrom_rename.txt --threads 4 -Oz -o cactus_renamed_chr1_2Mb.vcf.gz

# Subset to 10 samples + filter to SVs (50-50000bp)
echo "$REF_SAMPLES" | tr ' ' '\n' > all10.txt
echo "$TEST_SAMPLE" >> all10.txt
$BCF view -S all10.txt --force-samples cactus_renamed_chr1_2Mb.vcf.gz \
    | $BCF view -i 'abs(strlen(REF)-strlen(ALT))>=50 && abs(strlen(REF)-strlen(ALT))<=50000' \
        -Oz -o cactus_svs_10.vcf.gz --threads 4
$TABIX -p vcf cactus_svs_10.vcf.gz
echo "  cactus SVs (10 samples, Chr1:1-2Mb): $($BCF view -H cactus_svs_10.vcf.gz | wc -l) records"

echo ""
echo "[$(date)] STEP 2: subset GrENE-Net VCF to Chr1:1-2Mb for the same 10 samples"
# bgzip+index a small region of GrENE first if not done
if [ ! -f grene_chr1_2Mb.vcf.gz ]; then
    grep -E '^#|^1\t[12]?[0-9]{6,7}\t' $GRENE 2>/dev/null \
        | awk -F'\t' 'NR<=1000 || ($1=="1" && $2>=1 && $2<=2000000) || /^#/' \
        | head -200000 \
        > grene_chr1_2Mb.vcf
    # Better: use bcftools on the bgzipped file, but it's not bgzipped yet
    # Easier: just take the relevant region directly
    rm grene_chr1_2Mb.vcf
    awk -F'\t' '/^#/ {print; next} $1=="1" && $2>=1 && $2<=2000000 {print}' $GRENE \
        > grene_chr1_2Mb.vcf
    $BGZIP grene_chr1_2Mb.vcf
    $TABIX -p vcf grene_chr1_2Mb.vcf.gz
fi
echo "  GrENE Chr1:1-2Mb total records: $($BCF view -H grene_chr1_2Mb.vcf.gz | wc -l)"
$BCF view -S all10.txt --force-samples grene_chr1_2Mb.vcf.gz \
    -Oz -o grene_10.vcf.gz --threads 4
$TABIX -p vcf grene_10.vcf.gz

echo ""
echo "[$(date)] STEP 3: build reference panel (9 founders, SNPs+SVs)"
$BCF view -S <(echo "$REF_SAMPLES" | tr ' ' '\n') --force-samples cactus_svs_10.vcf.gz \
    -Oz -o ref_svs_9.vcf.gz --threads 4
$TABIX -p vcf ref_svs_9.vcf.gz
$BCF view -S <(echo "$REF_SAMPLES" | tr ' ' '\n') --force-samples grene_10.vcf.gz \
    -Oz -o ref_snps_9.vcf.gz --threads 4
$TABIX -p vcf ref_snps_9.vcf.gz

# Concat SNPs + SVs (sort by position), then drop records with any missing GT
# (Beagle's reference panel must be fully non-missing)
$BCF concat -a -d none ref_snps_9.vcf.gz ref_svs_9.vcf.gz --threads 4 \
    | $BCF sort \
    | $BCF view -e 'F_MISSING > 0' \
    -Oz -o ref_9.vcf.gz
$TABIX -p vcf ref_9.vcf.gz
echo "  ref_9 records (SNPs+SVs, no missing): $($BCF view -H ref_9.vcf.gz | wc -l)"

echo ""
echo "[$(date)] STEP 4: build target panel (1 founder, SNPs known + SVs masked)"
$BCF view -S <(echo "$TEST_SAMPLE") --force-samples cactus_svs_10.vcf.gz \
    -Oz -o target_svs_truth.vcf.gz --threads 4
$TABIX -p vcf target_svs_truth.vcf.gz
# Mask SV genotypes to ./. — this VCF is already SV-only, so we mask all
$BCF +setGT target_svs_truth.vcf.gz -Oz -o target_svs_masked.vcf.gz -- -t a -n .
$TABIX -p vcf target_svs_masked.vcf.gz
$BCF view -S <(echo "$TEST_SAMPLE") --force-samples grene_10.vcf.gz \
    -Oz -o target_snps_1.vcf.gz --threads 4
$TABIX -p vcf target_snps_1.vcf.gz

$BCF concat -a -d none target_snps_1.vcf.gz target_svs_masked.vcf.gz --threads 4 \
    | $BCF sort -Oz -o target_1.vcf.gz
$TABIX -p vcf target_1.vcf.gz
echo "  target_1 records (SNPs+SVs-masked): $($BCF view -H target_1.vcf.gz | wc -l)"

echo ""
echo "[$(date)] STEP 5: run Beagle imputation"
java -Xmx16g -jar $BEAGLE \
    ref=ref_9.vcf.gz \
    gt=target_1.vcf.gz \
    chrom=1 \
    out=imputed_1 \
    nthreads=4 \
    seed=42 \
    2>&1 | tail -30
$TABIX -p vcf imputed_1.vcf.gz

echo ""
echo "[$(date)] STEP 6: compare imputed vs truth"
$PYTHON <<PYEOF
import pysam
imp = pysam.VariantFile('imputed_1.vcf.gz')
truth = pysam.VariantFile('target_svs_truth.vcf.gz')
sample = '$TEST_SAMPLE'

truth_by_pos = {}
for rec in truth.fetch():
    sz = abs(len(rec.ref) - len(rec.alts[0]))
    if sz < 50 or sz > 50000: continue
    gt = rec.samples[sample]['GT']
    if gt is None or len(gt) == 0: continue
    truth_by_pos[(rec.chrom, rec.pos, rec.ref, rec.alts[0])] = gt[0]
print(f'  truth SVs in test region: {len(truth_by_pos):,}')

n_total, n_match, n_mismatch, n_missing_gt = 0, 0, 0, 0
for rec in imp.fetch():
    try:
        sz = abs(len(rec.ref) - len(rec.alts[0]))
    except Exception:
        continue
    if sz < 50 or sz > 50000: continue
    key = (rec.chrom, rec.pos, rec.ref, rec.alts[0])
    if key not in truth_by_pos: continue
    n_total += 1
    imp_gt = rec.samples[sample]['GT']
    if imp_gt is None or imp_gt[0] is None:
        n_missing_gt += 1; continue
    if imp_gt[0] == truth_by_pos[key]:
        n_match += 1
    else:
        n_mismatch += 1

print(f'  imputed records found: {n_total:,}')
print(f'  matches: {n_match:,}  mismatches: {n_mismatch:,}  missing: {n_missing_gt:,}')
if n_total > 0:
    print(f'  concordance: {n_match/n_total:.3f}')
    print(f'  (random baseline = 0.5 for biallelic)')
PYEOF

echo ""
echo "[$(date)] DONE — smoke test result printed above"
