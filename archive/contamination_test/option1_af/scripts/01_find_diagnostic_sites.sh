#!/bin/bash
#SBATCH --job-name=diag_sites
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-5
#SBATCH --requeue
#SBATCH --output=logs/diag_chr%a.out
#SBATCH --error=logs/diag_chr%a.err

# For one chrom of the 1141 panel, find SNPs where every 231 GrENE founder
# is hom-ref (no '1' in any 231 GT). Those are sites where the alt allele
# exists ONLY in the non-231 1001G accessions = diagnostic sites for
# contamination.
#
# Pure awk - avoids bcftools plugins (not installed) and tag-header issues.
#
# Outputs (per chrom):
#   sites/diag_chr<N>.tsv       chrom,pos,ref,alt,ac_non231,an_non231,af_non231,carriers
#   sites/diag_chr<N>.targets   chrom\tpos\tref,alt   (for bcftools mpileup -T)

mkdir -p logs
set -eo pipefail

CH="${1:-${SLURM_ARRAY_TASK_ID}}"

ROOT=/global/scratch/users/tbellg/hapfire_sv/contamination_test
PANEL="${ROOT}/vcf/1001G_80pilot_israel_regmap_overlapping_biallelic_chr${CH}_mac7.recode.vcf"
S231="/global/scratch/users/tbellg/hapfire_sv/data/vcf_samples_231.txt"
OUT_DIR="${ROOT}/option1_af/sites"

OUT_TSV="${OUT_DIR}/diag_chr${CH}.tsv"
OUT_TARGETS="${OUT_DIR}/diag_chr${CH}.targets"

echo "[$(date)] chr${CH}: scanning ${PANEL}"

awk -v S231="${S231}" -F'\t' '
  BEGIN {
    while ((getline line < S231) > 0) { in231[line] = 1 }
    close(S231)
    OFS = "\t"
  }
  /^##/ { next }
  /^#CHROM/ {
    for (i = 10; i <= NF; i++) {
      sample_at[i] = $i
      if ($i in in231) { idx231[i] = 1 } else { idxNon[i] = 1 }
    }
    print "chrom","pos","ref","alt","ac_non231","an_non231","af_non231","carriers"
    next
  }
  {
    # quick check: any alt allele among the 231? if yes, skip immediately.
    skip = 0
    for (i in idx231) {
      if ($i ~ /1/) { skip = 1; break }
    }
    if (skip) next

    # AC in 231 == 0. Compute AC/AN/AF in non-231 + carrier list.
    ac_non = 0; an_non = 0; carriers = ""
    for (i in idxNon) {
      gt = $i
      n = split(gt, a, "|")
      for (k = 1; k <= n; k++) {
        if (a[k] == "0" || a[k] == "1") { an_non++ }
        if (a[k] == "1") { ac_non++ }
      }
      if (gt ~ /1/) {
        smp = sample_at[i]
        carriers = (carriers == "") ? smp : carriers "," smp
      }
    }
    if (ac_non == 0) next   # at MAC>=7 panel filter, this should never trigger

    af = ac_non / an_non
    # VCF cols: $1=CHROM $2=POS $3=ID $4=REF $5=ALT
    print $1, $2, $4, $5, ac_non, an_non, af, carriers
  }
' "${PANEL}" > "${OUT_TSV}"

n=$(($(wc -l < "${OUT_TSV}") - 1))
echo "[$(date)] chr${CH}: ${n} diagnostic sites"

# bcftools mpileup -T expects tab-separated chrom\tpos\tref,alt
awk -F'\t' 'NR>1 {printf "%s\t%s\t%s,%s\n", $1, $2, $3, $4}' "${OUT_TSV}" > "${OUT_TARGETS}"
echo "[$(date)] chr${CH}: wrote ${OUT_TARGETS}"
