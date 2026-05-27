#!/bin/bash
#SBATCH --job-name=imp_loo
#SBATCH --partition=bse
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/loo_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/loo_%j.err

# Leave-one-out validation: for each of N founders selected from the 80 with
# both cactus + GrENE genotypes, mask the cactus SVs, impute, compare to truth.
#
# For prototype, pick 10 random founders. For full validation, all 80.

set -euo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/work
mkdir -p loo

BEAGLE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/beagle.jar
BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/tabix

# Pick 10 leave-out founders
N_LOO=${N_LOO:-10}
shuf --random-source=<(yes 42 | head -1000) overlap_samples.txt | head -$N_LOO > loo/loo_samples.txt
echo "Leave-out samples (n=$N_LOO):"
cat loo/loo_samples.txt

for SAMPLE in $(cat loo/loo_samples.txt); do
    echo ""
    echo "[$(date)] LOO for $SAMPLE"
    OUT_DIR=loo/$SAMPLE
    mkdir -p $OUT_DIR

    # 1. Reference: 79 founders (drop the leave-out)
    grep -v "^${SAMPLE}$" overlap_samples.txt > $OUT_DIR/ref_samples.txt
    $BCF view -S $OUT_DIR/ref_samples.txt --force-samples ref_80.vcf.gz \
        -Oz -o $OUT_DIR/ref_79.vcf.gz --threads 4
    $TABIX -p vcf $OUT_DIR/ref_79.vcf.gz

    # 2. Target: just the leave-out founder, with SVs masked to ./.
    #    Easy way: take ref_80, subset to just SAMPLE, then mask SV records
    $BCF view -S <(echo $SAMPLE) --force-samples ref_80.vcf.gz \
        -Oz -o $OUT_DIR/target_one.vcf.gz --threads 4
    # The target_one VCF (subset from ref_80) already has cactus SVs only.
    # Mask all genotypes to ./. since this whole VCF is SVs.
    $BCF +setGT $OUT_DIR/target_one.vcf.gz -Oz -o $OUT_DIR/target_masked.vcf.gz -- -t a -n .
    $TABIX -p vcf $OUT_DIR/target_masked.vcf.gz

    # 3. Beagle impute
    for chrom in 1 2 3 4 5; do
        java -Xmx48g -jar $BEAGLE \
            ref=$OUT_DIR/ref_79.vcf.gz \
            gt=$OUT_DIR/target_masked.vcf.gz \
            chrom=$chrom \
            out=$OUT_DIR/imp_chr$chrom \
            nthreads=8 \
            seed=42 \
            2>&1 | tail -5
        $TABIX -p vcf $OUT_DIR/imp_chr${chrom}.vcf.gz
    done
    $BCF concat $OUT_DIR/imp_chr*.vcf.gz -Oz -o $OUT_DIR/imp_all.vcf.gz --threads 4
    $TABIX -p vcf $OUT_DIR/imp_all.vcf.gz
    rm $OUT_DIR/imp_chr*.vcf.gz*

    # 4. Compare imputed vs truth for SV records
    /home/tbellagio/miniforge3/envs/hapfm/bin/python <<PYEOF
import pysam
imp = pysam.VariantFile('$OUT_DIR/imp_all.vcf.gz')
truth = pysam.VariantFile('cactus_svs_renamed.vcf.gz')

n_total, n_match, n_mismatch, n_missing = 0, 0, 0, 0
sv_records = {}
for rec in truth.fetch():
    if abs(len(rec.ref) - len(rec.alts[0])) < 50: continue  # skip non-SVs
    truth_gt = rec.samples['$SAMPLE']['GT']
    if truth_gt is None or len(truth_gt) == 0: continue
    sv_records[(rec.chrom, rec.pos)] = (rec.ref, rec.alts, truth_gt[0])

for rec in imp.fetch():
    if abs(len(rec.ref) - len(rec.alts[0])) < 50: continue
    key = (rec.chrom, rec.pos)
    if key not in sv_records: continue
    n_total += 1
    imp_gt = rec.samples['$SAMPLE']['GT']
    if imp_gt is None or imp_gt[0] is None:
        n_missing += 1; continue
    if imp_gt[0] == sv_records[key][2]:
        n_match += 1
    else:
        n_mismatch += 1

if n_total > 0:
    accuracy = n_match / n_total
    print(f'$SAMPLE: SV concordance = {n_match}/{n_total} = {accuracy:.3f}  '
          f'(mismatch={n_mismatch}, missing={n_missing})')
PYEOF
done
echo "[$(date)] DONE"
