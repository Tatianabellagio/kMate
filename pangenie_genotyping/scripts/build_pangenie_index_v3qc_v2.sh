#!/bin/bash
#SBATCH --job-name=pg_index_v3qc_v2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/pg_index_v3qc_v2_%j.out
#SBATCH --error=logs/pg_index_v3qc_v2_%j.err

# Build a fresh PanGenie-index from the v3qc-v2 VCF (the paper-quality, methodologically
# self-consistent version). Workflow:
#   1. Re-merge multi-allelics on v3qc_v2.haploid (bcftools norm -m +any)
#   2. Diploidize haploid GTs (0 → 0|0, 1 → 1|1, . → .|.) for PanGenie input
#   3. Run PanGenie-index → produces idx_Chr{1..5}_kmers.tsv.gz + Graph files
#
# Output prefix: pang_v3qc_v2_pangenie_index
# Mirrors the v3 build at /pangenie_genotyping/data/pang_135_pangenie_index_*.
mkdir -p logs
set -euo pipefail
eval "$(conda shell.bash hook)"
conda activate pangenie

BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping
SRC_HAP=$BASE/data/v3qc_v2/founders_231_v3qc_v2.haploid.vcf.gz
TMP=$BASE/data/v3qc_v2/idx_tmp
mkdir -p $TMP
OUT_PREFIX=$BASE/data/pang_v3qc_v2_pangenie_index

BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
PG_INDEX=/global/home/users/tbellg/miniforge3/envs/pangenie/bin/PanGenie-index

[ -s "$SRC_HAP" ] || { echo "ERROR: missing $SRC_HAP (haploidize_v3qc_v2 must run first)"; exit 1; }
[ -s "$REF" ] || { echo "ERROR: missing $REF"; exit 1; }

# ---- Step 1: re-merge multi-allelics so PanGenie sees them as bubbles ----
MERGED_MULTI=$TMP/founders_231_v3qc_v2.haploid.multiallelic.vcf.gz
if [ ! -s "$MERGED_MULTI" ]; then
    echo "[$(date)] Step 1: re-merge multi-allelics (norm -m +any)"
    $BCF norm -m +any $SRC_HAP --threads 4 -Oz -o $MERGED_MULTI 2> $TMP/multi_norm.log
    $BCF index -t $MERGED_MULTI 2>/dev/null || true
fi
N_HAP=$($BCF index -n $SRC_HAP)
N_MULTI=$($BCF index -n $MERGED_MULTI 2>/dev/null || $BCF view -H $MERGED_MULTI | wc -l)
echo "  records: $N_HAP biallelic → $N_MULTI multi-allelic-collapsed"

# ---- Step 2: diploidize for PanGenie (X → X|X) ----
PANG_VCF=$TMP/founders_231_v3qc_v2.diploid.vcf
if [ ! -s "$PANG_VCF" ]; then
    echo "[$(date)] Step 2: diploidize for PanGenie input"
    zcat $MERGED_MULTI | awk 'BEGIN{OFS="\t"} /^#/{print; next} {for(i=10;i<=NF;i++) $i=$i"|"$i; print}' > $PANG_VCF
fi
echo "  diploid VCF lines: $(wc -l < $PANG_VCF)"

# ---- Step 3: PanGenie-index ----
echo "[$(date)] Step 3: PanGenie-index"
$PG_INDEX -r $REF -v $PANG_VCF -o $OUT_PREFIX -t 8 -k 31

echo ""
echo "[$(date)] DONE"
ls -lh ${OUT_PREFIX}*
