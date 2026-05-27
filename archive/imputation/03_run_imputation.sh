#!/bin/bash
#SBATCH --job-name=imp_beagle
#SBATCH --partition=bse
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/beagle_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/beagle_%j.err

# Step 7 of the imputation plan: run Beagle to impute SVs into the 151 GrENE-only
# founders, using the 80-founder reference panel.

set -euo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/work

BEAGLE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/beagle.jar
TABIX=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/tabix

# Run per-chromosome (Beagle handles one chrom at a time most reliably)
for chrom in 1 2 3 4 5; do
    echo "[$(date)] Beagle on chrom $chrom"
    java -Xmx48g -jar $BEAGLE \
        ref=ref_80.vcf.gz \
        gt=target_151.vcf.gz \
        chrom=$chrom \
        out=imputed_151_chr${chrom} \
        nthreads=16 \
        seed=42 \
        2>&1 | tail -20
    $TABIX -p vcf imputed_151_chr${chrom}.vcf.gz
done

# Concatenate per-chrom outputs
echo "[$(date)] Concatenating chromosomes"
BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
$BCF concat \
    imputed_151_chr1.vcf.gz imputed_151_chr2.vcf.gz imputed_151_chr3.vcf.gz \
    imputed_151_chr4.vcf.gz imputed_151_chr5.vcf.gz \
    -Oz -o imputed_151.vcf.gz --threads 4
$TABIX -p vcf imputed_151.vcf.gz
rm imputed_151_chr*.vcf.gz*

echo "[$(date)] DONE"
echo "Imputed VCF: imputed_151.vcf.gz"
echo "  samples: $($BCF query -l imputed_151.vcf.gz | wc -l)"
echo "  records: $($BCF view -H imputed_151.vcf.gz | wc -l)"
