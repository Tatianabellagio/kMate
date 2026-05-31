#!/bin/bash
#SBATCH --job-name=kmer_pa_prod_v3
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --array=1-5
#SBATCH --output=logs/kmer_pa_prod_v3_%A_%a.out
#SBATCH --error=logs/kmer_pa_prod_v3_%A_%a.err
#
# PRODUCTION kmer_pa builder, v3qc_v3. One array task per chromosome (Chr1..Chr5).
# Generates the founder x k-mer matrix AND applies the production column filter
# IN THE SAME STEP (--filter-production: keep 2 <= ac <= F-1, dropping ac==0
# dead / ac==1 private / ac==F invariant). The matrix is written already-filtered
# -- no separate filt2 / filt2inv post-step. See ALGORITHM.md s2.1.
#
# Missing GTs: --treat-missing-as-n (N-on) is the production choice (./. -> REF
# would fabricate confident reference genotypes; see ALGORITHM.md s2.1 + M3).
#
# Condo QoS is sometimes CPU-capped (QOSGrpCpuLimit); if tasks pend, resubmit with
#   sbatch --qos=savio_lowprio --requeue scripts/build_kmer_pa_production_v3qc_v3.sh
# The skip-if-exists guard makes it safe to resubmit / requeue.
#
# NOTE: this writes the FILTERED production matrix. If you need the RAW matrix
# (all columns incl. invariant/dead) for a diagnostic -- e.g. the missingness
# test -- run build_kmer_pa.py WITHOUT --filter-production to a separate dir.
#
mkdir -p logs
set -euo pipefail
BASE=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python   # old 'hapfm' env is gone
$PY -c "import pysam, scipy, numpy" || { echo "ERROR: \$PY lacks pysam/scipy/numpy"; exit 1; }

CHR=Chr${SLURM_ARRAY_TASK_ID:-1}
KMERS=$BASE/panel/pangenie_genotyping/data/pang_135_pangenie_index_${CHR}_kmers.tsv.gz
VCF=$BASE/panel/pangenie_genotyping/data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
OUT_DIR=$BASE/data/kmer_pa_231_v3qc_v3_filt2inv          # production matrix dir
mkdir -p "$OUT_DIR"
OUT_PREFIX=$OUT_DIR/cn_${CHR}

[ -s "$KMERS" ] || { echo "ERROR: missing $KMERS"; exit 1; }
[ -s "$VCF" ]   || { echo "ERROR: missing $VCF"; exit 1; }
[ ! -s "${OUT_PREFIX}.kmer_pa.npz" ] || { echo "$CHR exists, skipping"; exit 0; }

echo "[$(date)] $CHR PRODUCTION kmer_pa (N-on + --filter-production)"
$PY -u $BASE/src/build_kmer_pa.py \
    --kmers "$KMERS" --vcf "$VCF" --ref "$REF" \
    --chrom "$CHR" --out "$OUT_PREFIX" \
    --treat-missing-as-n \
    --filter-production --min-ac 2 --invariant-margin 1
echo "[$(date)] $CHR DONE"
ls -lh ${OUT_PREFIX}.kmer_pa.npz ${OUT_PREFIX}.meta.npz
