#!/bin/bash
#SBATCH --job-name=p231_psd
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=logs/07h_psd_%j.out
#SBATCH --error=logs/07h_psd_%j.out

# =============================================================================
# Pool-size x depth accuracy sweep (Fig S13-style, hapFIRE paper), p231 arm.
# Sibling of p80/scripts/07j_run_kmate_poolsize_depth_p80.sh: production kMate
# (--unit chrom, per_founder normalize, filt2inv in-house index, kmer-weight
# uniform) across pool size N and depth COV. One run gives both the
# per-record AF tsv AND h_per_chrom.npz (no --h-only).
#
# Usage: sbatch 07h_run_kmate_poolsize_depth_p231.sh N COV SEED
#   N    = 2 | 5 | 20 | 50 | 150
#   COV  = 1 | 10
#   SEED = 42 | 43 | 44
# =============================================================================
mkdir -p logs
set -euo pipefail
N=${1:?Usage: N COV SEED}
COV=${2:?Usage: N COV SEED}
SEED=${3:?Usage: N COV SEED}
ROOT=/global/scratch/users/tbellg/kmate
CTRL=$ROOT/benchmarks/p231
PYTHON=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
DRIVER=$ROOT/src/per_sample_per_chrom.py

KMER=$ROOT/data/kmer_pa_231_arch3_filt2inv/kmer_pa
CN_VAR=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
CN_VAR_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz
CN_VAR_CALLED=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz

SUB=cov${COV}_n${N}_g0_s${SEED}_hotspots_p231_chr1
READS=$CTRL/sims/$SUB/reads
[ -s "$READS/r1.fq" ] || { echo "ERROR: missing $READS/r1.fq" >&2; exit 1; }
[ -s "${KMER}_Chr1.kmer_pa.npz" ] || { echo "ERROR: missing $KMER" >&2; exit 1; }

OUT_DIR=$CTRL/results/kmate_chrom_poolsize_depth/n${N}_cov${COV}_s${SEED}
mkdir -p $OUT_DIR
SAMPLE=p231_chrom_psd_n${N}_cov${COV}_s${SEED}
OUT_TSV=$OUT_DIR/${SAMPLE}.tsv

echo "[$(date)] p231 --unit chrom  N=$N  cov=$COV  seed=$SEED"
$PYTHON -u $DRIVER \
    --kmer-pa-prefix $KMER \
    --var-pa $CN_VAR --var-called $CN_VAR_CALLED --var-meta $CN_VAR_META \
    --reads $READS/r1.fq $READS/r2.fq \
    --sample $SAMPLE --out $OUT_TSV \
    --threads 8 --chroms Chr1 \
    --unit chrom --kmer-weight uniform

echo; echo "[$(date)] DONE — $OUT_TSV  (+ h_per_chrom.npz)"
