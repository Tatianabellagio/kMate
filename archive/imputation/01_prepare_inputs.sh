#!/bin/bash
#SBATCH --job-name=imp_prep
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/prep_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/prep_%j.err

# Step 1-4 of the imputation plan: prepare cactus-SV-only VCF with renamed
# samples and chromosomes, ready to merge with GrENE-Net SNPs.
#
# Outputs:
#   cactus_svs_renamed.vcf.gz  — 80 founders × cactus SVs only (SNPs dropped),
#                                samples renamed to 1001G IDs, chroms renamed to 1/2/3/4/5

set -euo pipefail
mkdir -p /carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/{work,logs}
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/work

# Use the BIALLELIC.NORM version (already produced by 56179). Beagle requires
# biallelic input — multi-allelic bubbles in the cactus filtered VCF would need
# splitting anyway.
CACTUS=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_test/raw_vcfbub_lv0_diploid.biallelic.norm.vcf.gz
BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/tabix

# 1. Build sample-rename map (Assembly_ID → 1001G_ID), only keeping samples that
#    are in GrENE-Net 231
PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python
$PYTHON <<'PYEOF' > sample_rename.txt
import pandas as pd
panel = pd.read_csv('/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/sv_panel_to_accession_id.tsv', sep='\t')
asm_to_1001g = dict(zip(panel.Assembly_ID.astype(str), panel.Accession_ID.astype(str)))
grene = set(open('/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/vcf_samples_231.txt').read().split())
import subprocess
res = subprocess.run('/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools query -l '
                     '/home/tbellagio/scratch/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.vcf.gz',
                     shell=True, capture_output=True, text=True)
for asm in res.stdout.strip().split():
    g = asm_to_1001g.get(asm)
    if g and g in grene:
        print(f'{asm}\t{g}')
PYEOF
echo "Sample rename map: $(wc -l < sample_rename.txt) entries"

# 2. Build keep-samples list (the 80 1001G IDs of cactus founders that ARE in GrENE)
awk '{print $2}' sample_rename.txt > keep_samples.txt
echo "Keep samples: $(wc -l < keep_samples.txt)"

# 3. Build chrom-rename map (Chr1 → 1, etc.)
echo -e "Chr1\t1\nChr2\t2\nChr3\t3\nChr4\t4\nChr5\t5" > chrom_rename.txt

# 4. Apply renames to cactus VCF, then subset to the 80 GrENE-overlap samples,
#    and keep only SVs (≥50bp REF/ALT length difference, INS/DEL only)
echo "[$(date)] Step 4: rename + subset + SV filter"
$BCF view $CACTUS -Oz -o tmp_renamed.vcf.gz \
    --threads 4
$BCF reheader --samples sample_rename.txt -o tmp_sample_renamed.vcf.gz tmp_renamed.vcf.gz
$BCF annotate --rename-chrs chrom_rename.txt tmp_sample_renamed.vcf.gz -Oz -o tmp_chrom_renamed.vcf.gz --threads 4
$TABIX -p vcf tmp_chrom_renamed.vcf.gz

# 5. Subset to 80 keep samples + filter to SVs (length diff ≥ 50bp).
#    Also exclude very-long SVs (>50kb) which slow Beagle and rarely impute well.
$BCF view -S keep_samples.txt --force-samples tmp_chrom_renamed.vcf.gz \
    | $BCF view -i 'abs(strlen(REF)-strlen(ALT))>=50 && abs(strlen(REF)-strlen(ALT))<=50000' \
        -Oz -o cactus_svs_renamed.vcf.gz --threads 4
$TABIX -p vcf cactus_svs_renamed.vcf.gz

echo "[$(date)] DONE"
echo "Final VCF: cactus_svs_renamed.vcf.gz"
echo "  samples: $($BCF query -l cactus_svs_renamed.vcf.gz | wc -l)"
echo "  records: $($BCF view -H cactus_svs_renamed.vcf.gz | wc -l)"
echo "  chromosomes: $($BCF view -H cactus_svs_renamed.vcf.gz | cut -f1 | sort -u | tr '\n' ' ')"
ls -la cactus_svs_renamed.vcf.gz*
rm -f tmp_renamed.vcf.gz tmp_sample_renamed.vcf.gz tmp_chrom_renamed.vcf.gz tmp_chrom_renamed.vcf.gz.tbi
