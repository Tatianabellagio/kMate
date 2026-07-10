#!/bin/bash
#SBATCH --job-name=vg_sv_map
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/vgsvmap_%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/vgsvmap_%x_%j.err
set -eo pipefail
# Graph-native SV pipeline part 1 (slow, once): vg giraffe -> GAM, then vg pack.
# vg call (fast, iterable) is run separately on the saved pack.
cd /global/scratch/users/tbellg/kmate
POOL=${1:?pool e.g. cov10_n80_g0_s42_hotspots_p80_chr1}
PANEL=${2:-p80}
T=${SLURM_CPUS_PER_TASK:-16}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vg
GDIR=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output
GBZ=$GDIR/pang_1001gplus_82acc.d2.gbz
DIST=$GDIR/pang_1001gplus_82acc.d2.dist
MIN=$GDIR/pang_1001gplus_82acc.d2.shortread.withzip.min
D=benchmarks/$PANEL/sims/$POOL
W=benchmarks/accuracy_vs_competitors/work
GAM=$W/${POOL}.gam
PACK=$W/${POOL}.pack

echo "[$(date)] giraffe -> GAM ($T threads) on $(hostname): $POOL"
vg giraffe -Z "$GBZ" -d "$DIST" -m "$MIN" --rescue-attempts 0 \
   -f "$D/reads/r1.fq" -f "$D/reads/r2.fq" -t "$T" > "$GAM"
echo "[$(date)] GAM done: $(ls -la $GAM | awk '{print $5}') bytes"
echo "[$(date)] vg pack (Q5)"
vg pack -x "$GBZ" -g "$GAM" -Q 5 -t "$T" -o "$PACK"
echo "[$(date)] DONE pack: $(ls -la $PACK | awk '{print $5}') bytes -> $PACK"
