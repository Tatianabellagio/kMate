#!/bin/bash
#SBATCH --job-name=imp_merge
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/merge_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/merge_%j.err

# Step 5-6 of the imputation plan:
#   - Build the reference panel: 80 founders × {GrENE SNPs ∪ cactus SVs}
#   - Build the target panel: 151 GrENE-only founders × GrENE SNPs only

set -euo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/work

BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/tabix
BGZIP=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bgzip

GRENE=/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_updatedVCF/greneNet_final_v1.1.recode.vcf

# 1. bgzip and index GrENE VCF if not done
if [ ! -f grene_all.vcf.gz ]; then
    echo "[$(date)] bgzip-ing GrENE-Net VCF (~3GB, ~10 min)"
    $BGZIP -c $GRENE > grene_all.vcf.gz
    $TABIX -p vcf grene_all.vcf.gz
fi

# 2. Define sample sets
#    overlap_samples.txt: 80 1001G IDs that have BOTH cactus + GrENE
awk '{print $2}' sample_rename.txt > overlap_samples.txt
#    not_in_cactus.txt: 151 GrENE 1001G IDs without cactus
$BCF query -l grene_all.vcf.gz > grene_all_samples.txt
grep -v -F -x -f overlap_samples.txt grene_all_samples.txt > not_in_cactus.txt || true
echo "Sample counts:"
echo "  overlap (cactus ∩ GrENE):    $(wc -l < overlap_samples.txt)"
echo "  GrENE-only (target):         $(wc -l < not_in_cactus.txt)"
echo "  total:                       $(wc -l < grene_all_samples.txt)"

# 3. Subset GrENE VCF to the 80 overlap samples (for the reference panel SNPs)
echo "[$(date)] Subsetting GrENE to 80 overlap samples"
$BCF view -S overlap_samples.txt --force-samples grene_all.vcf.gz \
    -Oz -o grene_80.vcf.gz --threads 4
$TABIX -p vcf grene_80.vcf.gz
echo "  records: $($BCF view -H grene_80.vcf.gz | wc -l)"

# 4. Subset GrENE VCF to the 151 GrENE-only samples (the target)
echo "[$(date)] Subsetting GrENE to 151 GrENE-only samples"
$BCF view -S not_in_cactus.txt --force-samples grene_all.vcf.gz \
    -Oz -o grene_151.vcf.gz --threads 4
$TABIX -p vcf grene_151.vcf.gz
echo "  records: $($BCF view -H grene_151.vcf.gz | wc -l)"

# 5. Concatenate (sort by position) GrENE SNPs (80) + cactus SVs (80) → reference panel
echo "[$(date)] Building reference panel (SNPs + SVs for the 80)"
# Concat GrENE SNPs + cactus SVs, sort, then filter records with any missing GT
# (Beagle's reference panel must be fully non-missing; smoke test confirmed this)
$BCF concat -a -d none --threads 4 \
    grene_80.vcf.gz cactus_svs_renamed.vcf.gz \
    | $BCF sort \
    | $BCF view -e 'F_MISSING > 0' \
    -Oz -o ref_80.vcf.gz
$TABIX -p vcf ref_80.vcf.gz
echo "  records: $($BCF view -H ref_80.vcf.gz | wc -l)"

# 6. Build target VCF: 151 GrENE-only founders, ALL variants but with SVs as ./.
#    Approach: subset cactus_svs_renamed to ZERO samples (header only), build a
#    new VCF body by adding 151 missing-GT columns. Then concat with grene_151.
echo "[$(date)] Building target panel (151 founders, SNPs only, SVs ./.)"
PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python
$PYTHON - <<'PYEOF'
import gzip, subprocess
BCF = "/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools"

# Read the 151 GrENE-only sample names
samples_151 = [s.strip() for s in open("not_in_cactus.txt")]
print(f"  building target with {len(samples_151)} samples")

# Open cactus_svs_renamed.vcf.gz, write a new VCF where each record's
# samples columns are replaced with 151 ./.
import pysam
src = pysam.VariantFile("cactus_svs_renamed.vcf.gz")
# Build new header with 151 sample names
hdr = pysam.VariantHeader()
for line in str(src.header).split("\n"):
    if line.startswith("##") and not line.startswith("##bcftools"):
        hdr.add_line(line)
for s in samples_151:
    hdr.add_sample(s)

with pysam.VariantFile("target_svs_missing.vcf.gz", "wz", header=hdr) as out:
    for rec in src.fetch():
        nr = out.new_record(contig=rec.contig, start=rec.start, stop=rec.stop,
                             id=rec.id, qual=rec.qual, filter=rec.filter,
                             alleles=(rec.ref, *rec.alts),
                             info=dict(rec.info))
        # All 151 samples → ./.
        for s in samples_151:
            nr.samples[s]["GT"] = (None, None)
        out.write(nr)
print(f"  wrote target_svs_missing.vcf.gz")
PYEOF
$TABIX -p vcf target_svs_missing.vcf.gz

# Concat with grene_151 → target_151 for imputation
$BCF concat -a -d none --threads 4 \
    grene_151.vcf.gz target_svs_missing.vcf.gz \
    | $BCF sort -Oz -o target_151.vcf.gz
$TABIX -p vcf target_151.vcf.gz

echo "[$(date)] DONE"
echo "Reference panel: ref_80.vcf.gz   ($($BCF query -l ref_80.vcf.gz | wc -l) samples)"
echo "Target panel:    target_151.vcf.gz ($($BCF query -l target_151.vcf.gz | wc -l) samples)"
ls -la ref_80.vcf.gz* target_151.vcf.gz*
