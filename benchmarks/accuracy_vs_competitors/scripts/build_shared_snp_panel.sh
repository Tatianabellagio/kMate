#!/bin/bash
# Build the ONE shared SNP panel that BOTH kMate and hapFIRE estimate AF for, so
# the accuracy comparison is fair (same sites, same founder genotypes). Source =
# the production arch3 231-founder panel VCF (merged_231_chr<N>_final.vcf.gz).
#
# Conventions (disclosed):
#  - biallelic SNPs only (ref_len==alt_len==1);
#  - DROP positions carrying >1 SNP record (split multiallelics) — hapFIRE/HARP
#    want one record per position;
#  - haploid founder GT -> phased homozygous diploid (inbred ecotypes): 0->0|0,
#    1->1|1; MISSING '.' -> 0|0 (impute REF) so hapFIRE sees a complete phased
#    panel (it errors on any unphased/missing GT). The MAR truth (built
#    separately) still counts only called founders, and we also report on the
#    high-call subset, so this impute choice does not silently bias the headline.
set -eo pipefail
ROOT=/global/scratch/users/tbellg/kmate
PANEL=${1:-p80}              # p80 (production target: all long-read founders) | p231
CHR=${2:-Chr1}
case "$PANEL" in
  p80)  SRC=$ROOT/benchmarks/p80/data/pangenome_${PANEL}_${CHR/Chr/chr}.vcf.gz ;;
  p231) SRC=$ROOT/panel/arch3/chr${CHR#Chr}/merged_231_${CHR/Chr/chr}_final.vcf.gz ;;
  *) echo "unknown panel: $PANEL"; exit 1 ;;
esac
WORK=$ROOT/benchmarks/accuracy_vs_competitors/work
OUT=$WORK/shared_snps_${PANEL}_${CHR}.vcf.gz
DUPS=$WORK/dup_pos_${PANEL}_${CHR}.txt

echo "[$(date)] build shared SNP panel: $SRC -> $OUT on $(hostname)"
test -f "$SRC" || { echo "missing source VCF: $SRC"; exit 1; }

# Pass 1 (cheap): positions of biallelic SNPs that occur >1x (split multiallelics).
bcftools view -H -v snps -m2 -M2 "$SRC" | cut -f2 \
 | sort -n | uniq -d > "$DUPS"
echo "  duplicate-POS to drop: $(wc -l < "$DUPS")"

# Pass 2: stream biallelic SNPs, skip dup positions, haploid GT -> phased homozygous.
bcftools view -v snps -m2 -M2 "$SRC" \
 | awk -v dupf="$DUPS" 'BEGIN{OFS="\t"; while((getline p < dupf)>0) dup[p]=1}
        /^#/ {print; next}
        { if(dup[$2]) next;
          for(j=10;j<=NF;j++){ g=$j;
            if(g=="0") $j="0|0"; else if(g=="1") $j="1|1"; else $j="0|0"; }
          print }' \
 | bcftools view -Oz -o "$OUT" -
bcftools index -t "$OUT"

echo "[$(date)] done. shared SNP records: $(bcftools index -n "$OUT")"
