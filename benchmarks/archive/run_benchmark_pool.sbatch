#!/bin/bash
#SBATCH --job-name=bench_pool
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=72G
#SBATCH --time=06:00:00
#SBATCH --output=benchmarks/scripts/logs/bench_%x_%j.out
#SBATCH --error=benchmarks/scripts/logs/bench_%x_%j.out
set -eo pipefail
# Unified per-pool kMate benchmark: count k-mers ONCE, run GLOBAL + BLOCK(dynld_K500),
# score both into a per-pool rows file (merged into benchmark_table.tsv afterward).
# Usage: sbatch run_benchmark_pool.sbatch <POOL> <PANEL>
cd /global/scratch/users/tbellg/kmate
mkdir -p benchmarks/scripts/logs benchmarks/benchmark_runs/tsv benchmarks/benchmark_runs/rows
POOL=${1:?pool}; PANEL=${2:?p80|p231}
T=${SLURM_CPUS_PER_TASK:-8}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate

# KP_GLOBAL: GLOBAL-mode kmer_pa. p231 (mixed panel) uses the REAL production
# filt2inv matrix. p80 (homogeneous, all-long-read control) uses NO filter --
# filt2/filt2inv + omega=1/m_b exist specifically to cancel the cactus/PG
# k-mer-count imbalance, which doesn't exist on this balanced panel
# (ALGORITHM.md M6; confirmed 2026-07-07).
# KP_BLOCK: WINDOW/BLOCK-mode kmer_pa, unchanged (not --local-only -> out of
# scope for the 2026-07-06/07 normalization-fix benchmark refresh).
case "$PANEL" in
  p231) KP_GLOBAL=data/kmer_pa_231_arch3_filt2inv/kmer_pa
        KP_BLOCK=data/kmer_pa_231_arch3_filt2inv/kmer_pa
        VP=panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
        VM=panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz ;;
  p80)  KP_GLOBAL=benchmarks/p80/data/kmer_pa_p80/kmer_pa
        KP_BLOCK=benchmarks/p80/data/kmer_pa_p80_filt2/kmer_pa
        VP=benchmarks/p80/data/var_pa_p80.var_pa.npz
        VM=benchmarks/p80/data/var_pa_p80.meta.npz ;;
  *) echo "bad panel $PANEL"; exit 1 ;;
esac
BLOCKS=analysis/grenenet_gea/blocks_mcf90/chr1_units_dynld_K500.tsv
SIMD=benchmarks/$PANEL/sims/$POOL
TRUTH=$SIMD/recomb_truth.tsv.gz
[ -f "$TRUTH" ] || TRUTH=$SIMD/recomb_truth_raw.tsv.gz   # p231 sims use _raw naming
R1=benchmarks/$PANEL/sims/$POOL/reads/r1.fq
R2=benchmarks/$PANEL/sims/$POOL/reads/r2.fq
test -f "$TRUTH" || { echo "missing truth: $TRUTH"; exit 1; }

OUTD=benchmarks/benchmark_runs/tsv
ROWS=benchmarks/benchmark_runs/rows/${POOL}.tsv
GLOBAL=$OUTD/${POOL}_global.tsv
BLOCK=$OUTD/${POOL}_block_dynldK500.tsv
JF=${TMPDIR:-/tmp}/${POOL}.jf          # node-local k-mer DB (count once)
rm -f "$ROWS"

echo "[$(date)] build-kmer-db (count once) on $(hostname): $POOL"
kmate build-kmer-db --reads "$R1" "$R2" --out "$JF" --threads "$T"

echo "[$(date)] kMate GLOBAL"
python -m kmate.per_sample_per_chrom --kmer-pa-prefix "$KP_GLOBAL" --var-pa "$VP" --var-meta "$VM" \
  --reads "$R1" "$R2" --kmer-db "$JF" --sample g_${POOL} --threads "$T" --chroms Chr1 \
  --block-mode global --kmer-weight uniform --out "$GLOBAL"

echo "[$(date)] kMate BLOCK (dynld_K500)"
python -m kmate.per_sample_per_chrom --kmer-pa-prefix "$KP_BLOCK" --var-pa "$VP" --var-meta "$VM" \
  --reads "$R1" "$R2" --kmer-db "$JF" --sample b_${POOL} --threads "$T" --chroms Chr1 \
  --block-mode window --kmer-weight inv_mb --blocks-tsv "$BLOCKS" --min-kmers-per-block 50 \
  --out "$BLOCK"

echo "[$(date)] score both -> $ROWS"
python3 benchmarks/build_benchmark_table.py --tool kMate --mode global \
  --pool "$POOL" --truth "$TRUTH" --est "$GLOBAL" --table "$ROWS"
python3 benchmarks/build_benchmark_table.py --tool kMate --mode block \
  --pool "$POOL" --truth "$TRUTH" --est "$BLOCK" --table "$ROWS"
rm -f "$JF"
echo "[$(date)] DONE $POOL"