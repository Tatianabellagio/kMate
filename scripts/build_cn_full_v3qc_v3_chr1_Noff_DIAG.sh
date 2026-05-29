#!/bin/bash
#SBATCH --job-name=cn_v3_Noff_diag
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=logs/cn_full_v3qc_v3_Noff_diag_%j.out
#SBATCH --error=logs/cn_full_v3qc_v3_Noff_diag_%j.err
#
# ============================ DIAGNOSTIC ONLY ============================
# This is NOT a production build. Production cn_full_v3qc_v3 is N-ON
# (--treat-missing-as-n; ./. -> N), which is the correct modeling choice
# (./.->REF would fabricate confident reference genotypes from low-quality
# no-calls). See ALGORITHM.md s2.1 + M3.
#
# This build is N-OFF (./.->REF, the build_kmer_cn.py default) for ONE purpose:
# to measure how much of the cactus/PG private-k-mer imbalance is caused by
# N-on dropping PG missing-GT k-mers vs. real biology. Compare its all-bubble
# private/side-only ratios against the N-on matrix and against the fully-called-
# bubble estimate from notebooks/MISSINGNESS_CAUSES_IMBALANCE.ipynb.
# Delete data/cn_full_231_v3qc_v3_Noff_diag/ once the question is answered.
# ========================================================================
mkdir -p logs
set -euo pipefail
BASE=/global/scratch/users/tbellg/kmate

# --- SET THIS to an env with pysam + scipy + numpy (the old 'hapfm' env is gone). ---
# e.g. PY=/global/home/users/tbellg/miniforge3/envs/<env>/bin/python
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
$PY -c "import pysam, scipy, numpy" || { echo "ERROR: \$PY lacks pysam/scipy/numpy; set PY to a suitable env"; exit 1; }

CHR=Chr1
KMERS=$BASE/panel/pangenie_genotyping/data/pang_135_pangenie_index_${CHR}_kmers.tsv.gz
VCF=$BASE/panel/pangenie_genotyping/data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
OUT_DIR=$BASE/data/cn_full_231_v3qc_v3_Noff_diag
mkdir -p $OUT_DIR
OUT_PREFIX=$OUT_DIR/cn_${CHR}

[ -s "$VCF" ] || { echo "ERROR: missing $VCF"; exit 1; }
[ ! -s "${OUT_PREFIX}.cn.npz" ] || { echo "exists, skipping"; exit 0; }

echo "[$(date)] $CHR cn_full v3qc_v3 N-OFF DIAGNOSTIC (./. -> REF; NO --treat-missing-as-n)"
$PY -u $BASE/src/build_kmer_cn.py \
    --kmers "$KMERS" --vcf "$VCF" --ref "$REF" \
    --chrom "$CHR" --out "$OUT_PREFIX"
echo "[$(date)] $CHR DONE"
ls -lh ${OUT_PREFIX}.cn.npz ${OUT_PREFIX}.meta.npz
