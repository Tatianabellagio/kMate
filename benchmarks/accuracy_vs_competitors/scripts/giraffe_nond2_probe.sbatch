#!/bin/bash
#SBATCH --job-name=gir_nd2p
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=0:40:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
set -uo pipefail
# Decide FAST whether the NON-d2 (unclipped) graph maps acceptably with node-local
# indexes + zipcode minimizer. 200k pairs only. If reads/s/thread is healthy
# (~1000s), the full end-to-end will sail. If it's pathological, pivot to d2.
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh; conda activate vg
T=${SLURM_CPUS_PER_TASK:-32}
GDIR=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output
PFX=pang_1001gplus_82acc
POOL=cov10_n50_g3_s42_hotspots_p80_chr1
R1=benchmarks/p80/sims/$POOL/reads/r1.fq; R2=benchmarks/p80/sims/$POOL/reads/r2.fq
LOC=${TMPDIR:-/tmp}/nd2probe; mkdir -p "$LOC"
echo "[$(date)] host=$(hostname) staging non-d2 indexes..."; t0=$SECONDS
cp $GDIR/$PFX.gbz $GDIR/$PFX.dist $GDIR/$PFX.shortread.withzip.min $GDIR/$PFX.shortread.zipcodes "$LOC/"
echo "[$(date)] staged in $((SECONDS-t0))s"
P=$LOC/$PFX
head -n 800000 $R1 > $LOC/r1.fq; head -n 800000 $R2 > $LOC/r2.fq   # 200k pairs
echo "[$(date)] giraffe non-d2 + zipcode (node-local), 200k pairs, --rescue-attempts 0"
/usr/bin/time -v vg giraffe -Z $P.gbz -d $P.dist -m $P.shortread.withzip.min -z $P.shortread.zipcodes \
   --rescue-attempts 0 -f $LOC/r1.fq -f $LOC/r2.fq -t $T -p --output-format GAF > $LOC/probe.gaf 2> $LOC/probe.err
echo "[$(date)] rc=$?"
grep -E "Mapping speed|Mapped .* reads|Elapsed \(wall|Maximum resident" $LOC/probe.err
echo "[$(date)] watchdog(>10s) reads: $(grep -c 'has been checked in for' $LOC/probe.err)"
echo "[$(date)] DONE"
