#!/usr/bin/env bash
# Extract per-record calling-confidence fields (F_MISSING, CONFLICT, MA) from the
# 231-founder Minigraph-Cactus panel VCFs, for non-SNP records only (ref_len!=1 or
# alt_len!=1) -- feeds the insertion-calling-artifact-vs-biology check (2026-07-03
# session, ancestral-polarization substitute: no outgroup genome/alignment exists in
# this repo, so we use the panel's OWN multi-assembly graph-conflict signal as a
# call-quality proxy instead of true ancestral state).
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
export PATH="/global/home/users/tbellg/miniforge3/envs/kmate/bin:$PATH"
OUT=analysis/grenenet_selection/sv_adaptive
mkdir -p "$OUT"
for c in 1 2 3 4 5; do
  echo "[chr$c] start $(date)"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/F_MISSING\t%INFO/CONFLICT\t%INFO/MA\n' \
    "panel/arch3/chr${c}/merged_231_chr${c}_final.vcf.gz" \
    | awk -F'\t' 'BEGIN{OFS="\t"} {
        rl=length($3); al=length($4);
        if (rl==1 && al==1) next;                 # SNP-only, drop
        conf = ($6=="."||$6=="") ? 0 : 1;
        fm = ($5=="."||$5=="") ? "nan" : $5;
        ma = ($7=="."||$7=="") ? "nan" : $7;
        print $1, $2, rl, al, fm, conf, ma
      }' \
    > "$OUT/vcf_callqual_chr${c}.tsv"
  echo "[chr$c] done $(date), $(wc -l < "$OUT/vcf_callqual_chr${c}.tsv") rows"
done
echo "ALL DONE $(date)"
