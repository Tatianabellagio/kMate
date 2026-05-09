#!/bin/bash
#SBATCH --job-name=mpileup_seedmix
#SBATCH --output=/home/tbellagio/scratch/hapfire_sv/contamination_test/option1_af/logs/mpileup_%x_%j.out
#SBATCH --error=/home/tbellagio/scratch/hapfire_sv/contamination_test/option1_af/logs/mpileup_%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# For one seedmix sample, run bcftools mpileup at the diagnostic sites
# (non-231-private SNPs, AC=0 in 231 GrENE founders) and write per-site
# allele depths to a TSV. AF = AD_alt / (AD_ref + AD_alt) at each site.
#
# Args:
#   $1 = sample N (1..8) -- maps to filtered_seeds-${N}.bam
#
# Output:
#   af/seedmix_S${N}.af.tsv       chrom,pos,ref,alt,ad_ref,ad_alt,dp,af

set -eo pipefail
source /home/tbellagio/miniforge3/etc/profile.d/conda.sh
conda activate pang
set -u

N="${1:?sample number 1-8 required}"

ROOT=/home/tbellagio/scratch/hapfire_sv/contamination_test/option1_af
BAM="/carnegie/nobackup/scratch/xwu/GrENE_net/seed_mix/filtered_bam/filtered_seeds-${N}.bam"
REF=/home/tbellagio/scratch/pang/ref_xing/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

# Targets file: chrom\tpos\tref,alt (built by 01_find_diagnostic_sites.sh).
# bcftools mpileup -T accepts this format directly.
TARGETS="${ROOT}/sites/diag_all.targets.gz"
[[ -f "${TARGETS}" ]] || { echo "missing targets: ${TARGETS}"; exit 2; }
n=$(zcat "${TARGETS}" | wc -l)
echo "[$(date)] targets: ${n} diagnostic sites"

OUT="${ROOT}/af/seedmix_S${N}.af.tsv"
mkdir -p "$(dirname "${OUT}")"

echo "[$(date)] mpileup S${N} at diagnostic sites..."

# -B: disable BAQ (faster, sufficient for AF estimation at known SNPs)
# -T: targets VCF (only emit records at these sites with these alts)
# -a FORMAT/AD,FORMAT/DP: emit allele depth and total depth
# --max-depth 1000: don't truncate (seedmix may have ~10-20x coverage genome-wide)
# -d 1000 isn't a flag here; we use --max-depth.
bcftools mpileup \
    -f "${REF}" \
    -T "${TARGETS}" \
    -a FORMAT/AD,FORMAT/DP \
    -B \
    --max-depth 1000 \
    --threads 4 \
    "${BAM}" 2>/dev/null \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\t[%DP]\n' \
  | awk 'BEGIN{FS=OFS="\t"; print "chrom","pos","ref","alt","ad_ref","ad_alt","dp","af"}
         {
           split($5, a, ",")
           ad_ref = (a[1] == "" || a[1] == ".") ? 0 : a[1]
           ad_alt = (a[2] == "" || a[2] == ".") ? 0 : a[2]
           dp = ($6 == "" || $6 == ".") ? 0 : $6
           af = (ad_ref + ad_alt > 0) ? ad_alt / (ad_ref + ad_alt) : 0
           print $1,$2,$3,$4,ad_ref,ad_alt,dp,af
         }' \
  > "${OUT}"

n=$(($(wc -l < "${OUT}") - 1))
echo "[$(date)] S${N}: wrote ${n} sites to ${OUT}"
