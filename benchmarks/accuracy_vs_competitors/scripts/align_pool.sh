#!/bin/bash
# Align a sim pool's reads to the TAIR10 linear reference for hapFIRE input.
# Replicates the p231 sim recipe (minimap2 -ax sr to TAIR10.chr.iupacN.fa) so the
# BAM contig naming (Chr1..Chr5) matches the shared SNP panel + reference fasta.
set -eo pipefail
ROOT=/global/scratch/users/tbellg/kmate
POOL=$1            # e.g. cov10_n231_g0_s42_hotspots_p80_chr1
PANEL=${2:-p80}
THREADS=${3:-8}
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
D=$ROOT/benchmarks/$PANEL/sims/$POOL
OUT=$ROOT/benchmarks/accuracy_vs_competitors/work/${POOL}.tair10.srt.bam

echo "[$(date)] align $POOL -> $OUT on $(hostname) ($THREADS threads)"
minimap2 -ax sr --MD --cs -Y --sam-hit-only -t "$THREADS" \
    -R '@RG\tID:illumina\tSM:bulk' \
    "$REF" "$D/reads/r1.fq" "$D/reads/r2.fq" \
 | samtools sort -@ 4 -o "$OUT" -
samtools index "$OUT"
echo "[$(date)] done. $(samtools idxstats "$OUT" | awk '{m+=$3} END{print m" mapped reads"}')"
