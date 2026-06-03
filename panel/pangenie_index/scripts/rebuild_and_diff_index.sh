#!/bin/bash
#SBATCH --job-name=idx_rebuild_diff
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=8:00:00
#SBATCH --output=logs/idx_rebuild_diff_%j.out
#SBATCH --error=logs/idx_rebuild_diff_%j.err
#
# LEVEL A (fresh): rebuild our in-house k-mer index FROM SCRATCH with
# build_kmers_tsv.py (no PanGenie), then diff against PanGenie's production index.
# Proves our builder reproduces PanGenie-index end-to-end this session.
#
# Same inputs as the original build (slurm_pang135_haploid.sh): pang_135 vcfbub
# VCF + non-iupacN TAIR10 ref, k=31, --haploid, jellyfish hash 3e9. Env: kmate
# (the removed 'pangenie' env is replaced; kmate has jellyfish + dna_jellyfish + pysam).
#
# Writes to a NEW dir (pang_135_haploid_rebuild) so it does NOT clobber the
# in-use ours_*.tsv.gz that the production K_pa builds consume.
#
mkdir -p logs
set -euo pipefail
BASE=/global/scratch/users/tbellg/kmate
WORK=$BASE/panel/pangenie_index
KENV=/global/home/users/tbellg/miniforge3/envs/kmate
export PATH=$KENV/bin:$PATH          # so build_kmers_tsv.py's `jellyfish` subprocess resolves to kmate
PY=$KENV/bin/python
$PY -c "import pysam, dna_jellyfish" || { echo "ERROR: kmate lacks pysam/dna_jellyfish"; exit 1; }
command -v jellyfish >/dev/null || { echo "ERROR: jellyfish not on PATH"; exit 1; }

VCF=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa
OUT_DIR=$WORK/pang_135_haploid_rebuild       # NEW dir — does not touch the in-use ours_*
OUT_PREFIX=$OUT_DIR/ours
PG_PREFIX=$BASE/panel/pangenie_genotyping/data/pang_135_pangenie_index
EXISTING_OURS=$WORK/pang_135_haploid/ours    # the index production currently consumes
mkdir -p "$OUT_DIR"

echo "[$(date)] === STAGE 1: rebuild our index from scratch (all 5 chroms) ==="
echo "  vcf=$VCF (haploid GT, raw vcfbub)"
echo "  ref=$REF (non-iupacN, matches original build)"
echo "  out=$OUT_PREFIX  (fresh; in-use ours_* untouched)"
/usr/bin/time -v $PY -u $WORK/scripts/build_kmers_tsv.py \
    --vcf "$VCF" --ref "$REF" --out "$OUT_PREFIX" \
    -k 31 --haploid --jellyfish-threads 8 --jellyfish-hash 3000000000

echo
echo "[$(date)] === STAGE 2a: diff FRESH-REBUILD vs PanGenie (Level A) ==="
$PY -u $WORK/scripts/diff_index_vs_pg.py \
    --pg "$PG_PREFIX" --ours "$OUT_PREFIX" --chroms Chr1,Chr2,Chr3,Chr4,Chr5 \
    | tee $BASE/results/index_levelA_freshrebuild_vs_pg.txt

echo
echo "[$(date)] === STAGE 2b: diff FRESH-REBUILD vs EXISTING ours_* (reproducibility) ==="
$PY -u $WORK/scripts/diff_index_vs_pg.py \
    --pg "$EXISTING_OURS" --ours "$OUT_PREFIX" --chroms Chr1,Chr2,Chr3,Chr4,Chr5 \
    | tee $BASE/results/index_levelA_freshrebuild_vs_existingours.txt

echo
echo "[$(date)] DONE — Level-A fresh rebuild + diffs in results/index_levelA_*.txt"
