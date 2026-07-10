#!/bin/bash
#SBATCH --job-name=gir_probe
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=0:30:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
set -uo pipefail
# Is giraffe slow because of the READS or the GRAPH? Map only the first 5k pairs on
# d2+zipcode with live -p progress (kept in the log) so we see the actual mapping RATE.
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh; conda activate vg
T=${SLURM_CPUS_PER_TASK:-16}
GDIR=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output
POOL=cov10_n50_g3_s42_hotspots_p80_chr1
SRC=benchmarks/p80/sims/$POOL/reads
TD=${TMPDIR:-/tmp}
head -20000 $SRC/r1.fq > $TD/r1.5k.fq   # 5000 reads
head -20000 $SRC/r2.fq > $TD/r2.5k.fq

echo "=== read length distribution (r1, first 5k) ==="
awk 'NR%4==2{print length($0)}' $TD/r1.5k.fq | sort -n | uniq -c | tail
echo "=== sample read names ==="
awk 'NR%4==1' $TD/r1.5k.fq | head -3

echo; echo "[$(date)] giraffe 5k pairs on d2+zip, T=$T, live progress:"
/usr/bin/time -v vg giraffe \
  -Z $GDIR/pang_1001gplus_82acc.d2.gbz \
  -d $GDIR/pang_1001gplus_82acc.d2.dist \
  -m $GDIR/pang_1001gplus_82acc.d2.shortread.withzip.min \
  -z $GDIR/pang_1001gplus_82acc.d2.shortread.zipcodes \
  -f $TD/r1.5k.fq -f $TD/r2.5k.fq -t $T -p --output-format GAF > $TD/probe.gaf
echo "[$(date)] mapped lines: $(wc -l < $TD/probe.gaf)"
echo "[$(date)] DONE"
