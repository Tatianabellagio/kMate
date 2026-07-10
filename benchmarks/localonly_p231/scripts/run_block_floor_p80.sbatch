#!/bin/bash
#SBATCH --job-name=blkflp80
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --output=benchmarks/localonly_p231/logs/%x_%A_%a.out
#SBATCH --error=benchmarks/localonly_p231/logs/%x_%A_%a.out
set -eo pipefail
# p80 per-block diagnostics: same floor=1 experiment as p231 but on the 80-founder
# (all-long-read) panel, in TWO k-mer-filter variants:
#   filt2      = private (ac=1) singletons DROPPED  (production-style)
#   unfiltered = private singletons KEPT             (the "delete the filter" test)
# Same 3 scenario analogs as the p231 run. w10kb unit (panel-agnostic window).
# Usage: sbatch --array=0-2 run_block_floor_p80.sbatch
cd /global/scratch/users/tbellg/kmate
mkdir -p benchmarks/localonly_p231/{floor_diag_p80,logs}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
export PYTHONPATH=/global/scratch/users/tbellg/kmate/src:$PYTHONPATH
T=${SLURM_CPUS_PER_TASK:-8}

POOLS=(
  cov10_n50_g0_s42_hotspots_p80_chr1
  cov10_n231_g1_s42_self97_hotspots_p80_chr1
  cov10_n50_g3_s42_hotspots_dom500nr_p80_chr1
)
POOL=${POOLS[$SLURM_ARRAY_TASK_ID]:?bad array index}

VP=benchmarks/p80/data/var_pa_p80.var_pa.npz
VM=benchmarks/p80/data/var_pa_p80.meta.npz
VC=benchmarks/p80/data/var_pa_p80.var_called.npz
SIMD=benchmarks/p80/sims/$POOL
R1=$SIMD/reads/r1.fq ; R2=$SIMD/reads/r2.fq
OUT=benchmarks/localonly_p231/floor_diag_p80
JF=${TMPDIR:-/tmp}/${POOL}.jf
DIAG=benchmarks/localonly_p231/scripts/block_floor_diag.py

echo "[$(date)] host=$(hostname) build-kmer-db: $POOL"
kmate build-kmer-db --reads "$R1" "$R2" --out "$JF" --threads "$T"

for VARIANT in filt2 unfiltered; do
  [ "$VARIANT" = filt2 ] && KP=benchmarks/p80/data/kmer_pa_p80_filt2/kmer_pa \
                         || KP=benchmarks/p80/data/kmer_pa_p80/kmer_pa
  echo "[$(date)] floor-diag w10kb variant=$VARIANT (KP=$KP)"
  python "$DIAG" --pool "$POOL" --unit w10kb --variant "$VARIANT" \
    --kmer-pa-prefix "$KP" --var-pa "$VP" --var-meta "$VM" --var-called "$VC" \
    --kmer-db "$JF" --chrom Chr1 --threads "$T" --outdir "$OUT"
done

rm -f "$JF"
echo "[$(date)] DONE $POOL"
