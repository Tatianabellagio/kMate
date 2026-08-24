#!/bin/bash
# Pangenome-native (founder-NAIVE) AF baseline: map a pool's reads to the p80
# cactus graph with vg giraffe, surject to TAIR10, then estimate per-SNP AF from
# allele depth at the shared SNP sites (AF = AD_alt / (AD_ref+AD_alt)). This is
# the "AF from coverage after graph mapping" comparator (cf. freqk paper's
# vg-giraffe -> pool-caller route); uses the SAME graph the p80 panel was
# deconstructed from, so variants/coords match the shared panel exactly.
set -eo pipefail
# Repo-root-relative: the kmate repo was relocated, so derive ROOT from this
# script's own location (benchmarks/accuracy_vs_competitors/scripts/) instead of
# hardcoding an absolute path that breaks on the next move.
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
POOL=${1:-cov10_n231_g0_s42_hotspots_p80_chr1}
PANEL=${2:-p80}
T=${3:-4}
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
WORK=$ROOT/benchmarks/accuracy_vs_competitors/work
SITES=$WORK/shared_snps_${PANEL}_Chr1.vcf.gz
D=$ROOT/benchmarks/$PANEL/sims/$POOL

# Graph selection by panel:
#  p80  -> external 82-accession cactus graph (ref paths named TAIR10#0#ChrN)
#  p231 -> vg autoindex construct graph from TAIR10 Chr1 + arch3 merged_231 SNP
#          VCF (build_giraffe_graph_p231.sbatch); ref path named plainly ChrN
ZIP=""
case "$PANEL" in
  p80)
    GDIR=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output
    GBZ=$GDIR/pang_1001gplus_82acc.d2.gbz
    DIST=$GDIR/pang_1001gplus_82acc.d2.dist
    MIN=$GDIR/pang_1001gplus_82acc.d2.shortread.withzip.min
    REFPATH_GREP='^TAIR10#0#Chr'
    ;;
  p231)
    IDXDIR=$WORK/autoidx_p231_snp_Chr1
    GBZ=$IDXDIR/graph.giraffe.gbz
    DIST=$IDXDIR/graph.dist
    MIN=$(ls "$IDXDIR"/graph*.min 2>/dev/null | head -1)
    ZIP=$(ls "$IDXDIR"/graph*.zipcodes 2>/dev/null | head -1 || true)
    REFPATH_GREP='^Chr'
    ;;
  *) echo "unknown panel: $PANEL" >&2; exit 1 ;;
esac
# --- Stage the giraffe index NODE-LOCAL to avoid Lustre mmap contention. ---
# Confirmed root cause of the p231 6h timeouts: giraffe memory-maps the index
# (esp. the ~526MB minimizer); one job maps this pool in 1:05, but 6 jobs reading
# the same shared-scratch index concurrently ALL exceeded 90min (even alone on a
# node -> not CPU/mem, it's shared-FS I/O). Copying to $TMPDIR (node-local disk)
# gives each job its own copy. cp is a one-time ~few-sec cost vs hours of contention.
STAGE=${TMPDIR:-/tmp}/vggraph.$$
mkdir -p "$STAGE"
trap 'rm -rf "$STAGE"' EXIT
echo "[$(date)] staging giraffe index -> $STAGE (node-local)"
for f in "$GBZ" "$DIST" "$MIN" ${ZIP:+"$ZIP"}; do cp "$f" "$STAGE/"; done
GBZ=$STAGE/$(basename "$GBZ"); DIST=$STAGE/$(basename "$DIST"); MIN=$STAGE/$(basename "$MIN")
[ -n "${ZIP:-}" ] && ZIP=$STAGE/$(basename "$ZIP")

PRESET=${4:-}                         # giraffe -b preset (e.g. fast); empty = default
PTAG=${PRESET:+_$PRESET}
BAM=$WORK/${POOL}.vg${PTAG}.tair10.srt.bam
OUT=$ROOT/benchmarks/accuracy_vs_competitors/results/vg_${POOL}${PTAG}_snp_frequency.txt

echo "[$(date)] vg giraffe map $POOL on $(hostname) ($T threads) preset='${PRESET:-default}'"
# ref-paths file: the TAIR10 linear reference paths to surject onto.
# MUST be per-job (node-local $STAGE): the content keys only on $PANEL, so a shared
# $WORK/${PANEL}_refpaths.txt is written+read by every concurrent pool at once ->
# one job truncates it while another reads -> giraffe "No sequence dictionary
# available" (killed 5/6 cov10 jobs in the first concurrent relaunch).
REFPATHS=$STAGE/${PANEL}_refpaths.txt
vg paths -L -x "$GBZ" | grep "$REFPATH_GREP" > "$REFPATHS"
echo "  surject onto $(wc -l < "$REFPATHS") ref paths (grep '$REFPATH_GREP')"
# 1) giraffe -> surjected BAM on TAIR10 reference paths (giraffe timed: it's the cost)
VGT=$ROOT/benchmarks/accuracy_vs_competitors/results/vg_${POOL}${PTAG}_time.txt
/usr/bin/time -v -o "$VGT" \
 vg giraffe -Z "$GBZ" -d "$DIST" -m "$MIN" ${ZIP:+-z "$ZIP"} ${PRESET:+-b "$PRESET"} \
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
