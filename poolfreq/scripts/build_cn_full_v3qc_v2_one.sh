#!/bin/bash
#SBATCH --job-name=cn_full_v3qc_v2
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/cn_full_v3qc_v2_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/cn_full_v3qc_v2_%A_%a.err

# Build cn_full_231_v3qc_v2 with --treat-missing-as-n on (./. → N, no REF default).
# Uses founders_231_v3qc_v2.haploid.vcf.gz as source.
set -euo pipefail
BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
CHR="Chr${SLURM_ARRAY_TASK_ID:-1}"

KMERS=$BASE/pangenie_genotyping/data/pang_135_pangenie_index_${CHR}_kmers.tsv.gz
VCF=$BASE/pangenie_genotyping/data/v3qc_v2/founders_231_v3qc_v2.haploid.vcf.gz
REF=/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
OUT_DIR=$BASE/poolfreq/data/cn_full_231_v3qc_v2
OUT_PREFIX=$OUT_DIR/cn_${CHR}
mkdir -p $OUT_DIR

[ -s "$VCF" ] || { echo "ERROR: missing $VCF"; exit 1; }
[ ! -s "${OUT_PREFIX}.cn.npz" ] || { echo "$CHR exists, skipping"; exit 0; }

echo "[$(date)] $CHR: build cn_full_v3qc_v2 (--treat-missing-as-n ON)"
$PY -u $BASE/poolfreq/src/build_kmer_cn.py \
    --kmers "$KMERS" --vcf "$VCF" --ref "$REF" \
    --chrom "$CHR" --out "$OUT_PREFIX" \
    --treat-missing-as-n

echo "[$(date)] $CHR: DONE"
ls -lh ${OUT_PREFIX}.cn.npz ${OUT_PREFIX}.meta.npz
$PY <<EOF
import scipy.sparse as sp
M = sp.load_npz("${OUT_PREFIX}.cn.npz")
print(f'  sanity: shape={M.shape} nnz={M.nnz:,} density={M.nnz/(M.shape[0]*M.shape[1])*100:.3f}%')
EOF
