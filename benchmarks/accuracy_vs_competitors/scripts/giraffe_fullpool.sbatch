#!/bin/bash
#SBATCH --job-name=gir_full
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=0:45:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
set -uo pipefail
# Full pool, ONE graph config, STANDALONE (no preceding load to thrash cache), live -p.
# Usage: sbatch giraffe_fullpool.sbatch <d2|nond2>
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh; conda activate vg
T=${SLURM_CPUS_PER_TASK:-32}
MODE=${1:?d2 or nond2}
GDIR=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output
POOL=cov10_n50_g3_s42_hotspots_p80_chr1
R1=benchmarks/p80/sims/$POOL/reads/r1.fq; R2=benchmarks/p80/sims/$POOL/reads/r2.fq
if [ "$MODE" = d2 ]; then P=$GDIR/pang_1001gplus_82acc.d2; M=$P.shortread.withzip.min; Z=$P.shortread.zipcodes
else P=$GDIR/pang_1001gplus_82acc; M=$P.shortread.withzip.min; Z=$P.shortread.zipcodes; fi
echo "[$(date)] MODE=$MODE host=$(hostname) T=$T pool=$POOL"
/usr/bin/time -v vg giraffe -Z $P.gbz -d $P.dist -m $M -z $Z \
  -f $R1 -f $R2 -t $T -p --output-format GAF > ${TMPDIR:-/tmp}/${MODE}_full.gaf
echo "[$(date)] mapped lines: $(wc -l < ${TMPDIR:-/tmp}/${MODE}_full.gaf)  MODE=$MODE DONE"
