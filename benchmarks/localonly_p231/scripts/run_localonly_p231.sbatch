#!/bin/bash
#SBATCH --job-name=lo_p231
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=05:00:00
#SBATCH --output=benchmarks/localonly_p231/logs/%x_%A_%a.out
#SBATCH --error=benchmarks/localonly_p231/logs/%x_%A_%a.out
set -eo pipefail
# Benchmark the global-free window mode (--local-only) on the p231 panel across
# every simulation setup (selfing + outcrossing, balanced + selection), seed s42.
# Per pool, count k-mers ONCE then run 5 modes and score each NaN-aware:
#   global                                       (chrom-wide h reference)
#   anchored / localonly  x  {dynldK500, w10kb}  (the fallback-vs-abstain contrast
#                                                 at the canonical LD block unit and
#                                                 at 10kb fixed windows)
# Usage: sbatch --array=0-10 run_localonly_p231.sbatch
cd /global/scratch/users/tbellg/kmate
mkdir -p benchmarks/localonly_p231/{work,results,logs}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
export PYTHONPATH=/global/scratch/users/tbellg/kmate/src:$PYTHONPATH
T=${SLURM_CPUS_PER_TASK:-8}

POOLS=(
  cov10_n231_g0_s42_hotspots_p231_chr1
  cov10_n231_g1_s42_hotspots_p231_chr1
  cov10_n231_g1_s42_self97_hotspots_p231_chr1
  cov10_n50_g0_s42_hotspots_p231_chr1
  cov10_n50_g1_s42_hotspots_p231_chr1
  cov10_n50_g1_s42_self97_hotspots_p231_chr1
  cov10_n50_g3_s42_hotspots_p231_chr1
  cov10_n50_g3_s42_self97_hotspots_p231_chr1
  cov10_n50_g3_s42_hotspots_dom500_p231_chr1
  cov10_n50_g3_s42_hotspots_dom500nr_p231_chr1
  cov10_n50_g3_s42_self97_hotspots_dom500_p231_chr1
)
POOL=${POOLS[$SLURM_ARRAY_TASK_ID]:?bad array index}

KP=data/kmer_pa_231_arch3_filt2inv/kmer_pa
VP=panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
VM=panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz
VC=panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz
BLOCKS=analysis/grenenet_gea/blocks_mcf90/chr1_units_dynld_K500.tsv
SIMD=benchmarks/p231/sims/$POOL
TRUTH=$SIMD/recomb_truth_raw.tsv.gz
R1=$SIMD/reads/r1.fq ; R2=$SIMD/reads/r2.fq
test -f "$TRUTH" || { echo "missing truth: $TRUTH"; exit 1; }

WORK=benchmarks/localonly_p231/work
ROWS=benchmarks/localonly_p231/results/rows_${POOL}.tsv
JF=${TMPDIR:-/tmp}/${POOL}.jf
SCORE=benchmarks/localonly_p231/scripts/score_localonly.py
rm -f "$ROWS"

echo "[$(date)] host=$(hostname) build-kmer-db (count once): $POOL"
kmate build-kmer-db --reads "$R1" "$R2" --out "$JF" --threads "$T"

# common kMate invocation; mode + kmer-weight + other flags appended by each call.
# normalize is left at the per_sample_per_chrom.py default (per_founder, the
# 2026-07-06 Kf_w fix) for every mode.
run() {  # <label> <out> <flags...>
  local label=$1 out=$2; shift 2
  echo "[$(date)] kMate $label"
  python -m kmate.per_sample_per_chrom \
    --kmer-pa-prefix "$KP" --var-pa "$VP" --var-meta "$VM" --var-called "$VC" \
    --reads "$R1" "$R2" --kmer-db "$JF" --sample "${label}_${POOL}" \
    --threads "$T" --chroms Chr1 --out "$out" "$@"
}

G=$WORK/${POOL}_global.tsv
AD=$WORK/${POOL}_anchored_dynldK500.tsv
LD=$WORK/${POOL}_localonly_dynldK500.tsv
AW=$WORK/${POOL}_anchored_w10kb.tsv
LW=$WORK/${POOL}_localonly_w10kb.tsv

# GLOBAL mode: drop omega=1/m_b (uniform) per PIPELINE_STATE.md Sec.0 (2026-07-06).
run global   "$G"  --block-mode global --kmer-weight uniform
# WINDOW mode: --kmer-weight inv_mb unchanged (production window-mode choice).
run anchD    "$AD" --block-mode window --blocks-tsv "$BLOCKS" --min-kmers-per-block 50 \
                   --global-anchor-weight 0.3 --hmm-smooth-passes 5 --kmer-weight inv_mb
run localD   "$LD" --block-mode window --blocks-tsv "$BLOCKS" --min-kmers-per-block 50 --local-only --kmer-weight inv_mb
run anchW    "$AW" --block-mode window --window-bp 10000 --min-kmers-per-block 50 \
                   --global-anchor-weight 0.3 --hmm-smooth-passes 5 --kmer-weight inv_mb
run localW   "$LW" --block-mode window --window-bp 10000 --min-kmers-per-block 50 --local-only --kmer-weight inv_mb

echo "[$(date)] score all 5 -> $ROWS"
score() { python3 "$SCORE" --pool "$POOL" --mode "$1" --unit "$2" \
            --truth "$TRUTH" --est "$3" --table "$ROWS"; }
score global    none      "$G"
score anchored  dynldK500 "$AD"
score localonly dynldK500 "$LD"
score anchored  w10kb     "$AW"
score localonly w10kb     "$LW"

rm -f "$JF"
echo "[$(date)] DONE $POOL"
