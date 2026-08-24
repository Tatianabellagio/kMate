#!/bin/bash
#SBATCH --job-name=seedmix_duprate
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/seedmix_duprate_%j.out
#SBATCH --error=logs/seedmix_duprate_%j.err
#
# Measure the COORDINATE-BASED duplicate rate of SEEDMIX S1 — the same kind of
# duplication xwu's canonical pipeline (bwa mem -> Picard MarkDuplicates -> remove)
# acted on. Answers "are SEEDMIX reads PCR-duplicate-free?" with a number.
#
# Aligns the FULL trimmed (pre-dedup) seeds-1 FASTQ with minimap2 -ax sr (full
# coverage matters: coincidental coordinate-duplicates scale with depth, so a
# subsample would understate the rate), then samtools markdup -s reports the
# duplicate fraction. Reference is our TAIR10.chr.fa (exact ref is immaterial for
# a dup-rate estimate). NOTHING in production is touched.
#
mkdir -p logs
set -euo pipefail
B=/global/home/users/tbellg/miniforge3/envs/kmate/bin
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa
RD=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix_trimmed
R1=$RD/seeds-1.1.fastq.gz
R2=$RD/seeds-1.2.fastq.gz
OUT=/global/scratch/users/tbellg/kmate/analysis/panel_qc/seedmix_duprate
mkdir -p "$OUT"
TMP=$OUT/tmp_s1
mkdir -p "$TMP"

echo "[$(date)] align FULL trimmed SEEDMIX seeds-1 (minimap2 -ax sr) + markdup"
echo "  R1=$R1"; echo "  R2=$R2"; echo "  REF=$REF"

# minimap2 -> collate -> fixmate (-m for markdup) -> position sort -> markdup -s
$B/minimap2 -ax sr -t 8 "$REF" "$R1" "$R2" 2>$OUT/minimap2.log \
  | $B/samtools collate -@ 8 -O -u -T "$TMP/collate" - \
  | $B/samtools fixmate -@ 8 -m -u - - \
  | $B/samtools sort -@ 8 -m 2G -u -T "$TMP/sort" - \
  | $B/samtools markdup -@ 8 -s -f "$OUT/s1_markdup_stats.txt" - "$OUT/s1.markdup.bam"

echo
echo "[$(date)] === DUPLICATE STATS (SEEDMIX S1, coordinate-based) ==="
cat "$OUT/s1_markdup_stats.txt"

# also print a clean percent line
$B/python - "$OUT/s1_markdup_stats.txt" <<'EOF' 2>/dev/null || true
import sys,re
t=open(sys.argv[1]).read()
def g(k):
    m=re.search(rf"{k}:\s*(\d+)",t); return int(m.group(1)) if m else None
dup=g("DUPLICATE TOTAL"); exam=g("READ") or g("EXAMINED")
if dup and exam: print(f"\n>>> SEEDMIX S1 duplicate rate = {100*dup/exam:.2f}%  ({dup:,}/{exam:,})")
EOF

rm -rf "$TMP"
echo "[$(date)] DONE — stats in $OUT/s1_markdup_stats.txt"
