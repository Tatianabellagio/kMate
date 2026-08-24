#!/bin/bash
# Split the 1141-panel VCF into per-chrom plain-text VCFs that hapFIRE.py can open.
# hapFIRE's vcf_processing.py opens with open(vcf, "r") and parses GT via regex
# ([0-9])\|([0-9]) -- so we must (a) decompress, (b) keep biallelic SNPs only,
# (c) drop any record with missing or unphased GTs, (d) ensure non-"." variant IDs.

set -eo pipefail

# Source: panel_1141_phased.vcf.gz built by build_merged_1141.sh.
# That is the merge of:
#   - Panel B (1135 fully-phased 1001G Beagle panel)
#   - 6 GrENE-Net extras (100001/100002/6939/9940/9977/9992) extracted from the
#     fully-phased 231-panel greneNet_final_v1.1.recode.vcf
# restricted to sites present in BOTH panels (so all 1141 samples have phased GTs).
SRC=/global/scratch/users/tbellg/hapfire_sv/contamination_test/vcf/panel_1141_phased.vcf.gz
OUT=/global/scratch/users/tbellg/hapfire_sv/contamination_test/vcf
CH="${1:?chrom number required}"

source /global/home/users/tbellg/miniforge3/etc/profile.d/conda.sh
conda activate pang
set -u

echo "[$(date)] splitting chr${CH} -> ${OUT}/panel_1141_chr${CH}.vcf"

# Panel B is fully phased + zero missing, so no awk repair is needed.
# Pipeline:
#   1. biallelic SNPs only
#   2. set variant IDs (hapFIRE rejects "."-id records)
#   3. (optional safety) awk drops any record where any GT does not match
#      hapFIRE's regex /([0-9])\|([0-9])/. Should be 0 dropped on this panel.
bcftools view "$SRC" -r "$CH" -m2 -M2 -v snps -Ou \
  | bcftools annotate --set-id +'%CHROM\_%POS' -Ov \
  | awk '
      BEGIN{FS=OFS="\t"; n_drop=0; n_keep=0}
      /^#/{print; next}
      {
        ok=1
        for (i=10; i<=NF; i++) {
          n = index($i, ":")
          gt = (n == 0) ? $i : substr($i, 1, n-1)
          if (gt !~ /^[0-9]\|[0-9]$/) { ok=0; break }
        }
        if (ok) { n_keep++; print } else { n_drop++ }
      }
      END {
        print "[awk] kept=" n_keep " dropped=" n_drop " (non-phased-biallelic GT)" > "/dev/stderr"
      }' \
  > "${OUT}/panel_1141_chr${CH}.vcf.tmp"

mv "${OUT}/panel_1141_chr${CH}.vcf.tmp" "${OUT}/panel_1141_chr${CH}.vcf"

n=$(grep -v -c "^#" "${OUT}/panel_1141_chr${CH}.vcf" || true)
echo "[$(date)] chr${CH} done: ${n} records"
