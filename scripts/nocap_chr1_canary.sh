#!/bin/bash
#SBATCH --job-name=nocap_chr1
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=logs/nocap_chr1_%j.out
#SBATCH --error=logs/nocap_chr1_%j.err
#
# NO-CAPS CHECK (Chr1 canary). Builds the in-house k-mer index for Chr1 with
# --no-caps (lifts per-allele biallelic=16, multiallelic=32, and per-bubble
# total max(nr_paths,301); overhang cap unchanged), then builds K_pa off it and
# diffs vs the capped production K_pa (data/kmer_pa_231_arch3_filt2inv).
#
# Same inputs/recipe as the production builds:
#   index: build_kmers_tsv.py  --haploid -k31  on the pang_135 vcfbub VCF (Chr1 subset)
#          REF = TAIR10.chr.fa  (non-iupacN, as the index build always uses)
#   K_pa : build_kmer_pa.py     merged_231_chr1_final.vcf.gz + TAIR10.chr.iupacN.fa
#          --treat-missing-as-n --filter-production --min-ac 2 --invariant-margin 1
# Writes to NEW dirs; touches nothing the production pipeline consumes.
#
mkdir -p logs
set -euo pipefail
BASE=/global/scratch/users/tbellg/kmate
KENV=/global/home/users/tbellg/miniforge3/envs/kmate
export PATH=$KENV/bin:$PATH          # build_kmers_tsv.py spawns `jellyfish` as subprocess
PY=$KENV/bin/python
BC=$KENV/bin/bcftools
$PY -c "import pysam, dna_jellyfish, scipy, numpy" || { echo "ERROR: kmate env incomplete"; exit 1; }
command -v jellyfish >/dev/null || { echo "ERROR: jellyfish not on PATH"; exit 1; }

CHR=Chr1
SRC_VCF=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.vcf.gz
REF_IDX=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa
REF_KPA=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
MERGED=$BASE/panel/arch3/chr1/merged_231_chr1_final.vcf.gz

IDX_DIR=$BASE/panel/pangenie_index/pang_135_haploid_nocap
KPA_DIR=$BASE/data/kmer_pa_231_arch3_nocap_filt2inv
CAPPED_KPA=$BASE/data/kmer_pa_231_arch3_filt2inv
mkdir -p "$IDX_DIR" "$KPA_DIR" "$BASE/results"

CHR_VCF=$IDX_DIR/pang_${CHR}.vcf.gz
IDX_PREFIX=$IDX_DIR/ours
KMERS=$IDX_DIR/ours_${CHR}_kmers.tsv.gz
KPA_PREFIX=$KPA_DIR/kmer_pa_${CHR}

echo "[$(date)] === STAGE 0: subset source VCF to $CHR ==="
if [ ! -s "$CHR_VCF" ]; then
  $BC view -r "$CHR" "$SRC_VCF" -Oz -o "$CHR_VCF"
  echo "  wrote $CHR_VCF ($($BC index -n "$SRC_VCF" 2>/dev/null || echo '?') src records; subsetting $CHR)"
else
  echo "  $CHR_VCF exists, reusing"
fi

echo
echo "[$(date)] === STAGE 1: build NO-CAP index for $CHR ==="
/usr/bin/time -v $PY -u $BASE/panel/pangenie_index/scripts/build_kmers_tsv.py \
    --vcf "$CHR_VCF" --ref "$REF_IDX" --out "$IDX_PREFIX" \
    -k 31 --haploid --no-caps \
    --jellyfish-threads 8 --jellyfish-hash 3000000000
echo "  index lines (bubbles): $(zcat "$KMERS" | tail -n +2 | wc -l)"

echo
echo "[$(date)] === STAGE 2: build K_pa from NO-CAP index ($CHR) ==="
/usr/bin/time -v $PY -u $BASE/src/build_kmer_pa.py \
    --kmers "$KMERS" --vcf "$MERGED" --ref "$REF_KPA" \
    --chrom "$CHR" --out "$KPA_PREFIX" \
    --treat-missing-as-n \
    --filter-production --min-ac 2 --invariant-margin 1
ls -lh ${KPA_PREFIX}.kmer_pa.npz ${KPA_PREFIX}.meta.npz

echo
echo "[$(date)] === STAGE 3: compare NO-CAP K_pa vs CAPPED production K_pa ==="
$PY -u $BASE/scripts/compare_kmer_pa_dirs.py \
    --a "$KPA_DIR"    --a-label nocap \
    --b "$CAPPED_KPA" --b-label capped \
    --chroms "$CHR" | tee $BASE/results/nocap_vs_capped_kmer_pa_chr1.txt

echo
echo "[$(date)] DONE — report in results/nocap_vs_capped_kmer_pa_chr1.txt"
