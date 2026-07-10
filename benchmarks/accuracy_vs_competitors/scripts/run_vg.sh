#!/bin/bash
# Pangenome-native (founder-NAIVE) AF baseline: map a pool's reads to the p80
# cactus graph with vg giraffe, surject to TAIR10, then estimate per-SNP AF from
# allele depth at the shared SNP sites (AF = AD_alt / (AD_ref+AD_alt)). This is
# the "AF from coverage after graph mapping" comparator (cf. freqk paper's
# vg-giraffe -> pool-caller route); uses the SAME graph the p80 panel was
# deconstructed from, so variants/coords match the shared panel exactly.
set -eo pipefail
ROOT=/global/scratch/users/tbellg/kmate
POOL=${1:-cov10_n231_g0_s42_hotspots_p80_chr1}
PANEL=${2:-p80}
T=${3:-4}
GDIR=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output
GBZ=$GDIR/pang_1001gplus_82acc.d2.gbz
DIST=$GDIR/pang_1001gplus_82acc.d2.dist
MIN=$GDIR/pang_1001gplus_82acc.d2.shortread.withzip.min
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
WORK=$ROOT/benchmarks/accuracy_vs_competitors/work
SITES=$WORK/shared_snps_${PANEL}_Chr1.vcf.gz
D=$ROOT/benchmarks/$PANEL/sims/$POOL
PRESET=${4:-}                         # giraffe -b preset (e.g. fast); empty = default
PTAG=${PRESET:+_$PRESET}
BAM=$WORK/${POOL}.vg${PTAG}.tair10.srt.bam
OUT=$ROOT/benchmarks/accuracy_vs_competitors/results/vg_${POOL}${PTAG}_snp_frequency.txt

echo "[$(date)] vg giraffe map $POOL on $(hostname) ($T threads) preset='${PRESET:-default}'"
# ref-paths file: the TAIR10 linear reference paths to surject onto
REFPATHS=$WORK/tair10_refpaths.txt
vg paths -L -x "$GBZ" | grep '^TAIR10#0#Chr' > "$REFPATHS"
echo "  surject onto $(wc -l < "$REFPATHS") TAIR10 ref paths"
# 1) giraffe -> surjected BAM on TAIR10 reference paths (giraffe timed: it's the cost)
VGT=$ROOT/benchmarks/accuracy_vs_competitors/results/vg_${POOL}${PTAG}_time.txt
/usr/bin/time -v -o "$VGT" \
 vg giraffe -Z "$GBZ" -d "$DIST" -m "$MIN" ${PRESET:+-b "$PRESET"} \
   ${VG_RESCUE:+--rescue-attempts "$VG_RESCUE"} \
   -f "$D/reads/r1.fq" -f "$D/reads/r2.fq" \
   -o BAM --ref-paths "$REFPATHS" -t "$T" \
 | samtools sort -@ "$T" -o "$BAM" -
samtools index "$BAM"
echo "[$(date)] mapped. BAM contigs:"; samtools view -H "$BAM" | grep -m3 '^@SQ'

# 2) reheader contigs TAIR10#0#ChrN -> ChrN if needed (match REF + sites VCF)
if samtools view -H "$BAM" | grep -q 'SN:TAIR10#0#'; then
  echo "[$(date)] reheadering TAIR10#0#ChrN -> ChrN"
  samtools view -H "$BAM" | sed -E 's/SN:TAIR10#0#/SN:/' > "$WORK/${POOL}.vg.hdr.sam"
  samtools reheader "$WORK/${POOL}.vg.hdr.sam" "$BAM" > "$BAM.re" && mv "$BAM.re" "$BAM"
  samtools index "$BAM"
fi

# 3) allele depth at shared SNP sites -> AF
echo "[$(date)] mpileup AD at shared SNP sites -> $OUT"
bcftools mpileup -f "$REF" -T "$SITES" -a FORMAT/AD -Ou "$BAM" 2>/dev/null \
 | bcftools query -f '%CHROM\t%POS\t[%AD]\n' \
 | awk 'BEGIN{OFS="\t"} {n=split($3,a,","); ref=a[1]; alt=(n>=2?a[2]:0); dp=ref+alt;
        print $1, $2, (dp>0 ? alt/dp : "nan")}' > "$OUT"
echo "[$(date)] done. $(wc -l < "$OUT") sites -> $OUT"
