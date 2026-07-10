#!/bin/bash
#SBATCH --job-name=gir_abtest
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=1:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
set -uo pipefail
# Decisive A/B: does feeding giraffe the ZIPCODE index fix the 30-130s/read wall?
# Same recombinant pool, giraffe-only, output discarded. Each config capped at 20min:
# a fixed config finishes in ~2min; a still-broken one gets killed so we learn fast.
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh; conda activate vg
T=${SLURM_CPUS_PER_TASK:-32}
GDIR=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output
POOL=cov10_n50_g3_s42_hotspots_p80_chr1
R1=benchmarks/p80/sims/$POOL/reads/r1.fq
R2=benchmarks/p80/sims/$POOL/reads/r2.fq
echo "[$(date)] host=$(hostname) cpus=$T pool=$POOL"
echo "[$(date)] reads: $(grep -c '^@' /dev/null 2>/dev/null; wc -l < $R1 | awk '{print $1/4}') pairs"

run () {
  local tag=$1; shift
  echo; echo "===== [$tag] ====="; echo "[$(date)] start: vg giraffe $*"
  /usr/bin/time -v timeout 1200 vg giraffe "$@" -f "$R1" -f "$R2" -t "$T" -p --output-format GAF > /dev/null 2> ${TMPDIR:-/tmp}/${tag}.time
  rc=$?
  echo "[$(date)] $tag rc=$rc (124=timeout/killed-at-20min)"
  grep -E "Elapsed \(wall|Maximum resident" ${TMPDIR:-/tmp}/${tag}.time
}

# A: NON-d2 graph + its zipcode index (cheap fix: snarls still match vg call, no scoring rework)
run A_nond2_zip -Z $GDIR/pang_1001gplus_82acc.gbz \
  -d $GDIR/pang_1001gplus_82acc.dist \
  -m $GDIR/pang_1001gplus_82acc.shortread.withzip.min \
  -z $GDIR/pang_1001gplus_82acc.shortread.zipcodes

# B: d2 CLIPPED graph + zipcodes (community-standard fast path; would need node remap for call)
run B_d2_zip -Z $GDIR/pang_1001gplus_82acc.d2.gbz \
  -d $GDIR/pang_1001gplus_82acc.d2.dist \
  -m $GDIR/pang_1001gplus_82acc.d2.shortread.withzip.min \
  -z $GDIR/pang_1001gplus_82acc.d2.shortread.zipcodes

echo; echo "[$(date)] DONE"
