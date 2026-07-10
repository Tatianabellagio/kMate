#!/bin/bash
#SBATCH --job-name=oh_time_chr1
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=4:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# ONE-OFF (2026-06-01): measure the overhang-stage cost in the in-house k-mer
# index build. Runs build_kmers_tsv.py on the SAME Chr1-subset VCF twice:
#   A) default            -> overhang SKIPPED (new behavior)
#   B) --emit-overhang    -> overhang COMPUTED (old behavior)
# Identical input, so the wall-time delta is purely the overhang stage.
# NOTE: Chr1-subset => "genome-wide" uniqueness is Chr1-only, so the ABSOLUTE
# k-mer set differs slightly from the production whole-genome build; that is
# irrelevant here (both runs use the same subset, the comparison is internal).
set -euo pipefail
KENV=/global/home/users/tbellg/miniforge3/envs/kmate
export PATH=$KENV/bin:$PATH
PY=$KENV/bin/python
WORK=/global/scratch/users/tbellg/kmate/panel/pangenie_index
VCF_ALL=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa
OUT=/global/scratch/users/tbellg/kmate/results/overhang_timing
mkdir -p "$OUT"

CHR1_VCF=$OUT/pang_135_chr1.vcf.gz
if [ ! -s "$CHR1_VCF" ]; then
  echo "[$(date)] subsetting VCF to Chr1 ..."
  bcftools view -r Chr1 "$VCF_ALL" -Oz -o "$CHR1_VCF"
fi
echo "[$(date)] Chr1 records: $(bcftools view -H "$CHR1_VCF" | wc -l)"

echo "[$(date)] === RUN A: default (overhang OFF) ==="
/usr/bin/time -v $PY -u "$WORK/scripts/build_kmers_tsv.py" \
    --vcf "$CHR1_VCF" --ref "$REF" --out "$OUT/noOH" \
    -k 31 --haploid --jellyfish-threads 8 --jellyfish-hash 3000000000 \
    2> "$OUT/timeA_noOH.txt"
echo "RUN A wall/mem:"; grep -E "Elapsed \(wall|Maximum resident" "$OUT/timeA_noOH.txt"

echo "[$(date)] === RUN B: --emit-overhang (overhang ON, old behavior) ==="
/usr/bin/time -v $PY -u "$WORK/scripts/build_kmers_tsv.py" \
    --vcf "$CHR1_VCF" --ref "$REF" --out "$OUT/withOH" \
    -k 31 --haploid --emit-overhang --jellyfish-threads 8 --jellyfish-hash 3000000000 \
    2> "$OUT/timeB_withOH.txt"
echo "RUN B wall/mem:"; grep -E "Elapsed \(wall|Maximum resident" "$OUT/timeB_withOH.txt"

echo "[$(date)] === SANITY: kept column (unique_kmers, cols 1-4) identical? ==="
if cmp -s <(zcat "$OUT/noOH_Chr1_kmers.tsv.gz" | cut -f1-4) \
          <(zcat "$OUT/withOH_Chr1_kmers.tsv.gz" | cut -f1-4); then
  echo "  cols 1-4 BYTE-IDENTICAL  (overhang change does not touch K_pa input)"
else
  echo "  cols 1-4 DIFFER  (UNEXPECTED — investigate)"
fi
echo "  col5 (overhang) distinct values in default run: $(zcat "$OUT/noOH_Chr1_kmers.tsv.gz" | tail -n+2 | cut -f5 | sort -u | tr '\n' ' ')"
echo "[$(date)] DONE"
