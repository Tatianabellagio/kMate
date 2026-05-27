#!/bin/bash
#SBATCH --job-name=merge_1141
#SBATCH --output=/home/tbellagio/scratch/hapfire_sv/contamination_test/logs/merge_1141_%j.out
#SBATCH --error=/home/tbellagio/scratch/hapfire_sv/contamination_test/logs/merge_1141_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

# Build a fully-phased 1141-ecotype panel by merging:
#   - Panel B (1135 ecotypes, fully phased Beagle output)
#   - the 6 GrENE-Net extras (100001/100002/6939/9940/9977/9992) extracted
#     from the 231-panel greneNet_final_v1.1 (also fully phased Beagle output)
#
# The two panels have 225 overlapping sample names, so we restrict the 231 panel
# to ONLY the 6 extras before merging (avoids bcftools merge's duplicate-name error).
#
# Post-merge filter F_MISSING==0 keeps only sites present in BOTH panels --
# at those sites all 1141 samples have phased GTs.

set -eo pipefail

source /home/tbellagio/miniforge3/etc/profile.d/conda.sh
conda activate pang
set -u

PANEL_B=/carnegie/nobackup/scratch/xwu/GrENE_net/vcf/1001_genomes_snps_missing0.8_merged_imputed_biallelic_named.vcf.gz
PANEL_231=/carnegie/nobackup/scratch/xwu/GrENE_net/greneNet_final_v1.1.recode.vcf
OUT_DIR=/home/tbellagio/scratch/hapfire_sv/contamination_test/vcf
EXTRAS=100001,100002,6939,9940,9977,9992

mkdir -p "${OUT_DIR}"

EXTRAS_VCF="${OUT_DIR}/extras_6.vcf.gz"
MERGED_VCF="${OUT_DIR}/panel_1141_phased.vcf.gz"

echo "[$(date)] step 1: extract 6 GrENE extras + bgzip + index"
bcftools view -s "${EXTRAS}" "${PANEL_231}" -Oz -o "${EXTRAS_VCF}"
bcftools index --tbi --threads 4 "${EXTRAS_VCF}"
n_extras=$(bcftools view -H "${EXTRAS_VCF}" | wc -l)
echo "[$(date)] extras VCF: 6 samples x ${n_extras} records"

echo "[$(date)] step 2: merge Panel B + 6 extras (-m all -> all 1141 in output)"
bcftools merge -m all --threads 4 "${PANEL_B}" "${EXTRAS_VCF}" -Oz -o "${MERGED_VCF}.unfiltered.vcf.gz"
bcftools index --tbi --threads 4 "${MERGED_VCF}.unfiltered.vcf.gz"
n_unf=$(bcftools view -H "${MERGED_VCF}.unfiltered.vcf.gz" | wc -l)
echo "[$(date)] merged unfiltered: 1141 samples x ${n_unf} records"

echo "[$(date)] step 3: keep only intersection sites (F_MISSING=0)"
bcftools view -e 'F_MISSING > 0' "${MERGED_VCF}.unfiltered.vcf.gz" -Oz -o "${MERGED_VCF}"
bcftools index --tbi --threads 4 "${MERGED_VCF}"
n_final=$(bcftools view -H "${MERGED_VCF}" | wc -l)
echo "[$(date)] merged filtered: 1141 samples x ${n_final} records (kept $(awk -v a=$n_final -v b=$n_unf 'BEGIN{printf "%.1f%%", 100*a/b}'))"

echo "[$(date)] sanity: phasing check on first 1k records of merged"
bcftools view -H "${MERGED_VCF}" | head -1000 | awk -F'\t' '
{ ph=0; uph=0; mis=0
  for(i=10;i<=NF;i++){split($i,a,":"); g=a[1]
    if(g~/\./) mis++; else if(g~/\|/) ph++; else if(g~/\//) uph++}
  tp+=ph; tu+=uph; tm+=mis; n++ }
END{ printf "  avg phased=%.0f, unphased=%.0f, missing=%.0f (per record, %d sampled)\n", tp/n, tu/n, tm/n, n }'

echo "[$(date)] sanity: chrom counts in merged"
bcftools index -s "${MERGED_VCF}"

echo "[$(date)] DONE: ${MERGED_VCF}"
ls -lh "${MERGED_VCF}"*
# clean up the unfiltered intermediate
rm -f "${MERGED_VCF}.unfiltered.vcf.gz" "${MERGED_VCF}.unfiltered.vcf.gz.tbi"
