#!/bin/bash
#SBATCH --job-name=kmate_seedmix_demo
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=2:00:00
#SBATCH --output=kmate_seedmix_demo_%j.out
#SBATCH --error=kmate_seedmix_demo_%j.err
# -----------------------------------------------------------------------------
# kMate real-data example: estimate allele frequencies in a REAL pool.
#
# Runs kMate on a real GrENE-Net SEEDMIX pool (a sequenced mixture of the 231
# Arabidopsis founders) against the production 231-founder arch3 panel, on Chr1,
# in global mode (SEEDMIX is a founder/F0 seed pool -> one mixture per chrom).
#
# This is the "Tier 2" example: unlike `kmate selftest` (a tiny bundled fixture),
# it uses the LARGE real panel matrices + a multi-GB real read pool that live on
# the cluster (not in the git repo). Treat the paths below as the template to
# adapt to your own panel + reads.
#
# Memory: the dense Chr1 kmer_pa for 231 founders needs ~80 GB (hence --mem=80G).
# Submit with sbatch; do NOT run on a login node.
#
#   sbatch run_seedmix.sh            # full sample S1
#   SAMPLE=S2 sbatch run_seedmix.sh  # a different SEEDMIX sample
# -----------------------------------------------------------------------------
set -euo pipefail
hostname

ROOT=/global/scratch/users/tbellg/kmate
S=${SAMPLE:-S1}                                   # SEEDMIX sample S1..S8
OUT_DIR=${OUT_DIR:-$ROOT/examples/real_seedmix/results}
mkdir -p "$OUT_DIR"

# Activate the kmate env (must have the `kmate` command installed: pip install -e .)
source ~/miniforge3/etc/profile.d/conda.sh
conda activate kmate

# --- inputs (real, on-cluster) -----------------------------------------------
R1=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/${S}-1.1_P.fq.gz
R2=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/${S}-1.2_P.fq.gz
KMER_PA=$ROOT/data/kmer_pa_231_arch3_filt2inv/kmer_pa          # prefix; _Chr1.* appended
VAR=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
VAR_CALLED=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz
VAR_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz
for f in "$R1" "$R2" "$VAR" "$VAR_META" "${KMER_PA}_Chr1.kmer_pa.npz"; do
    [ -s "$f" ] || { echo "ERROR: missing input $f" >&2; exit 1; }
done

# --- run: estimate h + project to per-record AF (production recipe) -----------
# --kmer-weight inv_mb is the production per-bubble de-replication weight.
kmate run \
    --kmer-pa-prefix "$KMER_PA" \
    --var-pa "$VAR" --var-called "$VAR_CALLED" --var-meta "$VAR_META" \
    --reads "$R1" "$R2" \
    --sample "SEEDMIX_${S}" --out "$OUT_DIR/SEEDMIX_${S}_chr1.tsv" \
    --threads 4 --chroms Chr1 --block-mode global --kmer-weight inv_mb

echo "DONE -> $OUT_DIR/SEEDMIX_${S}_chr1.tsv (+ .h_per_chrom.npz)"
echo "Validate the recovered founder mixture against hapFIRE truth with:"
echo "  python $ROOT/scripts/validate_seedmix_vs_hapfire.py   # see README.md"
