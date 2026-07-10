#!/bin/bash
#SBATCH --job-name=gir_loc
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=0:50:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
set -uo pipefail
# Localize the giraffe slowdown. Three timed runs, node-LOCAL indexes (Lustre removed as a variable):
#   1) 50k pairs  paired   -> baseline per-read time
#   2) 200k pairs paired   -> if per-read time GROWS vs 50k => superlinear (memory/rescue accumulation)
#   3) 200k pairs r1-only  -> if MUCH faster than (2) => pair-rescue is the cost
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh; conda activate vg
T=${SLURM_CPUS_PER_TASK:-32}
GDIR=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output
SRC=$GDIR/pang_1001gplus_82acc.d2
POOL=cov10_n50_g3_s42_hotspots_p80_chr1
R1=benchmarks/p80/sims/$POOL/reads/r1.fq
R2=benchmarks/p80/sims/$POOL/reads/r2.fq

LOC=${TMPDIR:-/tmp}/giridx
mkdir -p $LOC
echo "[$(date)] host=$(hostname) T=$T  copying indexes to node-local $LOC ..."
t0=$SECONDS
cp $SRC.gbz $SRC.dist $SRC.shortread.withzip.min $SRC.shortread.zipcodes $LOC/
echo "[$(date)] copy done in $((SECONDS-t0))s ; df of local:"; df -hT $LOC | tail -1
P=$LOC/pang_1001gplus_82acc.d2

# node-local read subsets (reads are randomly ordered => head = random sample)
head -n 200000 $R1 > $LOC/r1_50k.fq;  head -n 200000 $R2 > $LOC/r2_50k.fq      # 50k pairs
head -n 800000 $R1 > $LOC/r1_200k.fq; head -n 800000 $R2 > $LOC/r2_200k.fq     # 200k pairs

run () {  # name  outgaf  extra-args...
  local name=$1; shift
  local gaf=$1; shift
  echo "================ RUN: $name ================"
  echo "[$(date)] start $name"
  /usr/bin/time -v vg giraffe -Z $P.gbz -d $P.dist -m $P.shortread.withzip.min -z $P.shortread.zipcodes \
      -t $T -p "$@" > $gaf 2> $LOC/$name.err
  local rc=$?
  local wall=$(grep -oE "Elapsed \(wall clock\) time.*" $LOC/$name.err | grep -oE "[0-9:.]+$")
  local maxrss=$(grep -oE "Maximum resident set size \(kbytes\): [0-9]+" $LOC/$name.err | grep -oE "[0-9]+$")
  local lines=$(wc -l < $gaf)
  echo "[$(date)] $name rc=$rc  wall=$wall  maxRSS=$(awk -v k=$maxrss 'BEGIN{printf "%.1fG", k/1048576}')  mapped_records=$lines"
  echo "[$(date)] $name reads/sec(total) and per-read times: see below"
  grep -E "reads per second|Mapped .* reads|too-large rescue|Watchdog" $LOC/$name.err | head -5
  # count reads that tripped the 10s watchdog
  echo "[$(date)] $name watchdog(>10s) reads: $(grep -c 'has been checked in for' $LOC/$name.err)"
}

run paired50k   $LOC/p50.gaf   -f $LOC/r1_50k.fq  -f $LOC/r2_50k.fq
run paired200k  $LOC/p200.gaf  -f $LOC/r1_200k.fq -f $LOC/r2_200k.fq
run single200k  $LOC/s200.gaf  -f $LOC/r1_200k.fq

echo "================ SUMMARY ================"
echo "Compare wall/per-read across the three. Linear+fast => env/centromere-tail; superlinear => mem/rescue accumulation; single<<paired => pair-rescue is the cost."
echo "[$(date)] DONE"
