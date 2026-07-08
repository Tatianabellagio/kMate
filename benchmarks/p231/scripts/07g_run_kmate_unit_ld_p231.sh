#!/bin/bash
#SBATCH --job-name=p231_ld
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=logs/07g_ld_%j.out
#SBATCH --error=logs/07g_ld_%j.out

# =============================================================================
# benchmarks/p231 — NEW corrected-GLOBAL production estimator (2026-07-07):
#   * PRODUCTION in-house-index panel  data/kmer_pa_231_arch3_filt2inv  (not the
#     benchmark-local pang_135 rebuild) → true end-to-end production benchmark.
#   * --unit ld --ld-r2 0.1  (r²-LD CompleteLDPartition blocks; the corrected
#     "global" estimator) + haploblock collapse (eps=0) + per_founder + uniform.
#   * LD blocks precomputed from the SAME arch3 var_pa (data/ld_blocks_r2_0.10_Chr1.tsv).
# Supersedes 07f (--block-mode global on the benchmark-local panel).
#
# Usage: sbatch 07g_run_kmate_unit_ld_p231.sh REGIME CNVAR
#   REGIME = n50_g0 n231_g0 n50_g1 n231_g1 n50_g3 n50_g3_dom500 (+ *_self97 …)
#   CNVAR  = raw | atomized
# =============================================================================
mkdir -p logs
set -euo pipefail
REGIME=${1:?Usage: REGIME CNVAR}
CNVAR=${2:?Usage: REGIME CNVAR}
[[ "$CNVAR" == "atomized" || "$CNVAR" == "raw" ]] || { echo "ERROR: CNVAR must be raw|atomized" >&2; exit 1; }

ROOT=/global/scratch/users/tbellg/kmate
CTRL=$ROOT/benchmarks/p231
PYTHON=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
DRIVER=$ROOT/src/per_sample_per_chrom.py
COV=10; SEED=42

case "$REGIME" in
    n50_g0)        SUBDIR="cov10_n50_g0_s42_hotspots_p231_chr1" ;;
    n231_g0)       SUBDIR="cov10_n231_g0_s42_hotspots_p231_chr1" ;;
    n50_g1)        SUBDIR="cov10_n50_g1_s42_hotspots_p231_chr1" ;;
    n231_g1)       SUBDIR="cov10_n231_g1_s42_hotspots_p231_chr1" ;;
    n50_g3)        SUBDIR="cov10_n50_g3_s42_hotspots_p231_chr1" ;;
    n50_g1_self97) SUBDIR="cov10_n50_g1_s42_self97_hotspots_p231_chr1" ;;
    n231_g1_self97) SUBDIR="cov10_n231_g1_s42_self97_hotspots_p231_chr1" ;;
    n50_g3_self97) SUBDIR="cov10_n50_g3_s42_self97_hotspots_p231_chr1" ;;
    n50_g3_dom500) SUBDIR="cov10_n50_g3_s42_hotspots_dom500_p231_chr1" ;;
    n50_g3_dom500nr) SUBDIR="cov10_n50_g3_s42_hotspots_dom500nr_p231_chr1" ;;
    n50_g3_dom500_self97) SUBDIR="cov10_n50_g3_s42_self97_hotspots_dom500_p231_chr1" ;;
    *) echo "ERROR: unknown REGIME '$REGIME'" >&2; exit 1 ;;
esac

WORK=$CTRL/sims/$SUBDIR
READS_DIR=$WORK/reads
[ -s "$READS_DIR/r1.fq" ] || { echo "ERROR: missing $READS_DIR/r1.fq — run 06 first" >&2; exit 1; }

# PRODUCTION in-house index panel (Chr1)
CN_KMER_PREFIX=$ROOT/data/kmer_pa_231_arch3_filt2inv/kmer_pa
[ -s "${CN_KMER_PREFIX}_Chr1.kmer_pa.npz" ] || { echo "ERROR: missing production kmer_pa $CN_KMER_PREFIX" >&2; exit 1; }

LD_BLOCKS=$CTRL/data/ld_blocks_r2_0.10_Chr1.tsv
[ -s "$LD_BLOCKS" ] || { echo "ERROR: missing $LD_BLOCKS" >&2; exit 1; }

if [ "$CNVAR" = "atomized" ]; then
    CN_VAR=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1_atomized.var_pa.npz
    CN_VAR_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1_atomized.meta.npz
    CN_VAR_CALLED=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1_atomized.var_called.npz
else
    CN_VAR=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
    CN_VAR_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz
    CN_VAR_CALLED=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz
fi

OUT_DIR=$CTRL/results/kmate_ldr01_${CNVAR}/${REGIME}
mkdir -p $OUT_DIR
SAMPLE=p231_ldr01_${CNVAR}_${REGIME}_cov${COV}_s${SEED}
OUT_TSV=$OUT_DIR/${SAMPLE}.tsv

echo "[$(date)] kMate --unit ld (r2=0.1) GLOBAL-prod  regime=$REGIME  var_pa=$CNVAR"
echo "  kmer_pa:   $CN_KMER_PREFIX  (production in-house index)"
echo "  ld_blocks: $LD_BLOCKS"
echo "  var_pa:    $CN_VAR"
echo "  out:       $OUT_TSV"

$PYTHON -u $DRIVER \
    --kmer-pa-prefix $CN_KMER_PREFIX \
    --var-pa $CN_VAR \
    --var-called $CN_VAR_CALLED \
    --var-meta $CN_VAR_META \
    --reads $READS_DIR/r1.fq $READS_DIR/r2.fq \
    --sample $SAMPLE \
    --out $OUT_TSV \
    --threads 8 --chroms Chr1 \
    --unit ld --ld-blocks $LD_BLOCKS --kmer-weight uniform

echo; echo "[$(date)] DONE — $OUT_TSV"; ls -lh $OUT_TSV
