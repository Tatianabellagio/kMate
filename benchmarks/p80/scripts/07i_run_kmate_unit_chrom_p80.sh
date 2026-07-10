#!/bin/bash
#SBATCH --job-name=p80_chrom
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=logs/07i_chrom_%j.out
#SBATCH --error=logs/07i_chrom_%j.out

# =============================================================================
# p80 control under --unit chrom (settled selfing/inbred PRODUCTION estimator),
# in-house-index panel. Sibling of 07h (--unit ld); fills the gap left when p80
# never got a real chrom run — plot_scenario_grids_ldr01.py was silently
# plotting empty placeholder cells for KMATE_UNITTAG=chrom because
# results/kmate_chrom_p80_<ARM>/ never existed.
#
# Two k-mer arms (ARM): raw | filt2inv — confirms the private-k-mer drop is
# ~neutral on the homogeneous panel.
#
# Usage: sbatch 07i_run_kmate_unit_chrom_p80.sh REGIME ARM
#   REGIME = n50_g0 n200_g0 n231_g0 n50_g1 n200_g1 n231_g1 n50_g3 n50_g3_dom500 n80_g0
#            (+ *_self97 selfing variants)
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
    n50_g1_self97)  SUB="cov10_n50_g1_s42_self97_hotspots_p80_chr1" ;;
    n231_g1_self97) SUB="cov10_n231_g1_s42_self97_hotspots_p80_chr1" ;;
    n50_g3_self97)  SUB="cov10_n50_g3_s42_self97_hotspots_p80_chr1" ;;
    n50_g3_dom500_self97) SUB="cov10_n50_g3_s42_self97_hotspots_dom500_p80_chr1" ;;
    *) echo "ERROR: unknown REGIME '$REGIME'" >&2; exit 1 ;;
esac

READS=$CTRL/sims/$SUB/reads
[ -s "$READS/r1.fq" ] || { echo "ERROR: missing $READS/r1.fq" >&2; exit 1; }
[ -s "${KMER}_Chr1.kmer_pa.npz" ] || { echo "ERROR: missing $KMER — run 03d build first" >&2; exit 1; }

VAR=$CTRL/data/var_pa_p80
OUT_DIR=$CTRL/results/kmate_chrom_p80_${ARM}/${REGIME}
mkdir -p $OUT_DIR
SAMPLE=p80_chrom_${ARM}_${REGIME}_cov10_s42
OUT_TSV=$OUT_DIR/${SAMPLE}.tsv

echo "[$(date)] p80 --unit chrom  regime=$REGIME  arm=$ARM"
echo "  kmer_pa: $KMER"
$PYTHON -u $DRIVER \
    --kmer-pa-prefix $KMER \
    --var-pa ${VAR}.var_pa.npz --var-called ${VAR}.var_called.npz --var-meta ${VAR}.meta.npz \
    --reads $READS/r1.fq $READS/r2.fq \
    --sample $SAMPLE --out $OUT_TSV \
    --threads 8 --chroms Chr1 \
    --unit chrom --kmer-weight uniform
echo; echo "[$(date)] DONE — $OUT_TSV"; ls -lh $OUT_TSV
