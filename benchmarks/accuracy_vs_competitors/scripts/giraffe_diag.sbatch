#!/bin/bash
#SBATCH --job-name=gir_diag
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=1:40:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
set -uo pipefail
# Diagnose the per-read explosion on d2+zip. ARG: default | bounded
#   default  -> stock params, persist GAF + capture watchdog so we can measure %done + slow-read coords
#   bounded  -> cap DP / hits to bound worst-case per-read time (candidate fix)
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh; conda activate vg
T=${SLURM_CPUS_PER_TASK:-32}
MODE=${1:?default|bounded}
GDIR=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output
P=$GDIR/pang_1001gplus_82acc.d2
POOL=cov10_n50_g3_s42_hotspots_p80_chr1
R1=benchmarks/p80/sims/$POOL/reads/r1.fq; R2=benchmarks/p80/sims/$POOL/reads/r2.fq
TOT=$(( $(wc -l < $R1) / 4 ))
OUT=benchmarks/accuracy_vs_competitors/work/diag_${MODE}.gaf   # persistent (not TMPDIR)
ERR=benchmarks/accuracy_vs_competitors/logs/diag_${MODE}_watchdog.err
echo "[$(date)] MODE=$MODE host=$(hostname) T=$T total_pairs=$TOT"

EXTRA=""
[ "$MODE" = bounded ] && EXTRA="--max-dp-cells 2000000 --hard-hit-cap 200 --max-extensions 400 --max-alignments 4"
echo "[$(date)] EXTRA='$EXTRA'"

/usr/bin/time -v vg giraffe -Z $P.gbz -d $P.dist -m $P.shortread.withzip.min -z $P.shortread.zipcodes \
  $EXTRA -f $R1 -f $R2 -t $T -p --output-format GAF > $OUT 2> $ERR
rc=$?
echo "[$(date)] rc=$rc"
echo "[$(date)] mapped lines: $(wc -l < $OUT) / $TOT pairs ($(awk -v a=$(wc -l < $OUT) -v b=$TOT 'BEGIN{printf "%.1f", 100*a/(2*b)}')% of reads)"
grep -E "Mapping speed|Mapped .* reads across|Elapsed \(wall" $ERR 2>/dev/null
echo "[$(date)] slow-read coords (watchdog 'checked in', first field after origin):"
grep -oE "checked in for [0-9]+ seconds processing: c1h1fh_[^|]+\|[0-9]+" $ERR 2>/dev/null | grep -oE "[0-9]+$" | sort -n | awk '{print int($1/1000000)"Mb"}' | sort | uniq -c | sort -rn | head -20
echo "[$(date)] total watchdog-slow reads: $(grep -c 'has been checked in for' $ERR 2>/dev/null)"
echo "[$(date)] DONE MODE=$MODE"
