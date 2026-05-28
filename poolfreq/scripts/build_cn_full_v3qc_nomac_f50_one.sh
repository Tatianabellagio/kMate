#!/bin/bash
#SBATCH --job-name=cn_full_nomac_f50
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=logs/cn_full_nomac_f50_%A_%a.out
#SBATCH --error=logs/cn_full_nomac_f50_%A_%a.err

# Build cn_full from no-MAC + F_MISSING<=0.5 VCF for one chrom.
mkdir -p logs
set -euo pipefail
BASE=/global/scratch/users/tbellg/hapfire_sv
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CHR="Chr${SLURM_ARRAY_TASK_ID:-1}"

KMERS=$BASE/pangenie_genotyping/data/pang_135_pangenie_index_${CHR}_kmers.tsv.gz
VCF=$BASE/pangenie_genotyping/data/v3qc/founders_231_v3qc_noMAC_F50.haploid.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
OUT_DIR=$BASE/poolfreq/data/cn_full_231_v3qc_noMAC_F50
OUT_PREFIX=$OUT_DIR/cn_${CHR}
mkdir -p $OUT_DIR

[ -s "$VCF" ] || { echo "ERROR: missing $VCF"; exit 1; }
[ ! -s "${OUT_PREFIX}.cn.npz" ] || { echo "$CHR exists, skipping"; exit 0; }

echo "[$(date)] $CHR: build cn_full_v3qc_noMAC_F50"
$PY -u $BASE/poolfreq/src/build_kmer_cn.py \
    --kmers "$KMERS" --vcf "$VCF" --ref "$REF" \
    --chrom "$CHR" --out "$OUT_PREFIX"

echo "[$(date)] $CHR: DONE"
ls -lh ${OUT_PREFIX}.cn.npz ${OUT_PREFIX}.meta.npz
$PY <<EOF
import scipy.sparse as sp
M = sp.load_npz("${OUT_PREFIX}.cn.npz")
print(f'  sanity: shape={M.shape} nnz={M.nnz:,} density={M.nnz/(M.shape[0]*M.shape[1])*100:.3f}%')
EOF
