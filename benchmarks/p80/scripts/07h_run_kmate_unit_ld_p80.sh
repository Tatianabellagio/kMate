#!/bin/bash
#SBATCH --job-name=p80_ld
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=logs/07h_ld_%j.out
#SBATCH --error=logs/07h_ld_%j.out

# =============================================================================
# p80 control under the NEW --unit ld estimator, in-house-index panel.
# Two k-mer arms (ARM): raw | filt2inv — confirms the private-k-mer drop is
# ~neutral on the homogeneous panel. --unit ld --ld-r2 0.1 (LD blocks from
# var_pa_p80) + per_founder + uniform.
#
# Usage: sbatch 07h_run_kmate_unit_ld_p80.sh REGIME ARM
#   REGIME = n50_g0 n200_g0 n231_g0 n50_g1 n200_g1 n231_g1 n50_g3 n50_g3_dom500 n80_g0
#   ARM    = raw | filt2inv
# =============================================================================
mkdir -p logs
set -euo pipefail
REGIME=${1:?Usage: REGIME ARM}
ARM=${2:?Usage: REGIME ARM}
ROOT=/global/scratch/users/tbellg/kmate
CTRL=$ROOT/benchmarks/p80
PYTHON=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
DRIVER=$ROOT/src/per_sample_per_chrom.py

case "$ARM" in
    raw)      KMER=$CTRL/data/kmer_pa_p80_ours/kmer_pa ;;
    filt2inv) KMER=$CTRL/data/kmer_pa_p80_ours_filt2inv/kmer_pa ;;
    *) echo "ERROR: ARM must be raw|filt2inv" >&2; exit 1 ;;
esac
case "$REGIME" in
    n50_g0)  SUB="cov10_n50_g0_s42_hotspots_p80_chr1" ;;
    n200_g0) SUB="cov10_n200_g0_s42_hotspots_p80_chr1" ;;
    n231_g0) SUB="cov10_n231_g0_s42_hotspots_p80_chr1" ;;
    n80_g0)  SUB="cov10_n80_g0_s42_hotspots_p80_chr1" ;;
    n50_g1)  SUB="cov10_n50_g1_s42_hotspots_p80_chr1" ;;
    n200_g1) SUB="cov10_n200_g1_s42_hotspots_p80_chr1" ;;
    n231_g1) SUB="cov10_n231_g1_s42_hotspots_p80_chr1" ;;
    n50_g3)  SUB="cov10_n50_g3_s42_hotspots_p80_chr1" ;;
    n50_g3_dom500) SUB="cov10_n50_g3_s42_hotspots_dom500_p80_chr1" ;;
    *) echo "ERROR: unknown REGIME '$REGIME'" >&2; exit 1 ;;
esac

READS=$CTRL/sims/$SUB/reads
[ -s "$READS/r1.fq" ] || { echo "ERROR: missing $READS/r1.fq" >&2; exit 1; }
[ -s "${KMER}_Chr1.kmer_pa.npz" ] || { echo "ERROR: missing $KMER — run 03d build first" >&2; exit 1; }

VAR=$CTRL/data/var_pa_p80
OUT_DIR=$CTRL/results/kmate_ldr01_p80_${ARM}/${REGIME}
mkdir -p $OUT_DIR
SAMPLE=p80_ldr01_${ARM}_${REGIME}_cov10_s42
OUT_TSV=$OUT_DIR/${SAMPLE}.tsv

echo "[$(date)] p80 --unit ld  regime=$REGIME  arm=$ARM"
echo "  kmer_pa: $KMER"
$PYTHON -u $DRIVER \
    --kmer-pa-prefix $KMER \
    --var-pa ${VAR}.var_pa.npz --var-called ${VAR}.var_called.npz --var-meta ${VAR}.meta.npz \
    --reads $READS/r1.fq $READS/r2.fq \
    --sample $SAMPLE --out $OUT_TSV \
    --threads 8 --chroms Chr1 \
    --unit ld --ld-r2 0.1 --kmer-weight uniform
echo; echo "[$(date)] DONE — $OUT_TSV"; ls -lh $OUT_TSV
