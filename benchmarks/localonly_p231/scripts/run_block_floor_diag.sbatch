#!/bin/bash
#SBATCH --job-name=blkfloor
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --output=benchmarks/localonly_p231/logs/%x_%A_%a.out
#SBATCH --error=benchmarks/localonly_p231/logs/%x_%A_%a.out
set -eo pipefail
# Per-block identifiability diagnostics to justify/replace the --min-kmers-per-block
# floor. 3 representative pools span the founder-count (50 vs 231 = identifiability)
# and regime (g0 baseline / selfing / recombinant+selection) axes.
# Each task: build kmer-db once, then floor=1 no-smooth localonly fit for BOTH units.
# Usage: sbatch --array=0-2 run_block_floor_diag.sbatch
cd /global/scratch/users/tbellg/kmate
mkdir -p benchmarks/localonly_p231/{floor_diag,logs}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
export PYTHONPATH=/global/scratch/users/tbellg/kmate/src:$PYTHONPATH
T=${SLURM_CPUS_PER_TASK:-8}

POOLS=(
  cov10_n50_g0_s42_hotspots_p231_chr1            # outcross g0 baseline (global is "correct")
  cov10_n231_g1_s42_self97_hotspots_p231_chr1    # selfing, full 231 panel (hardest identifiability)
  cov10_n50_g3_s42_hotspots_dom500nr_p231_chr1   # recombinant + selection (local mode's use case)
)
POOL=${POOLS[$SLURM_ARRAY_TASK_ID]:?bad array index}

KP=data/kmer_pa_231_arch3_filt2inv/kmer_pa
VP=panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
VM=panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz
VC=panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz
BLOCKS=analysis/grenenet_gea/blocks_mcf90/chr1_units_dynld_K500.tsv
SIMD=benchmarks/p231/sims/$POOL
R1=$SIMD/reads/r1.fq ; R2=$SIMD/reads/r2.fq
OUT=benchmarks/localonly_p231/floor_diag
JF=${TMPDIR:-/tmp}/${POOL}.jf
DIAG=benchmarks/localonly_p231/scripts/block_floor_diag.py

echo "[$(date)] host=$(hostname) build-kmer-db: $POOL"
kmate build-kmer-db --reads "$R1" "$R2" --out "$JF" --threads "$T"

for UNIT in w10kb dynldK500; do
  echo "[$(date)] floor-diag $UNIT"
  EXTRA=""; [ "$UNIT" = dynldK500 ] && EXTRA="--blocks-tsv $BLOCKS"
  python "$DIAG" --pool "$POOL" --unit "$UNIT" \
    --kmer-pa-prefix "$KP" --var-pa "$VP" --var-meta "$VM" --var-called "$VC" \
    --kmer-db "$JF" --chrom Chr1 --threads "$T" --outdir "$OUT" $EXTRA
done

rm -f "$JF"
echo "[$(date)] DONE $POOL"
