#!/bin/bash
# ROADMAP_GLOBAL_REFRESH Phase 5 — re-score the 4-tool table with the REFRESHED
# kMate arm (--unit chrom, --normalize per_founder, --kmer-weight uniform, the
# in-house-index p80 panel kmer_pa_p80_ours_filt2inv).
#
# Sibling of score_all_competitors.sh, which stays as the pre-refresh record.
# Differences, all deliberate:
#
#  1. kMate est comes from benchmarks/p80/results/kmate_chrom_p80_filt2inv/<regime>/,
#     NOT benchmark_runs/tsv/. The latter was last written 2026-07-07 08:44 — after
#     per_founder (a8ba02d, 07-06 22:06) but BEFORE the full-panel Kf_w fix
#     (9669be7, 07-07 10:22), haploblock collapse (b4d6ce0, 16:15) and the --unit
#     unification (a9bf1c0, 18:12). Those files are stale for this purpose.
#
#  2. NO kMate-block arm. block_dynldK500 is window mode, which ROADMAP_GLOBAL_REFRESH
#     puts explicitly out of scope ("do not touch the 10 kb / star2 window recipe").
#     Emitting a refreshed chrom row beside a pre-refresh block row in one table would
#     silently mix two estimator versions, so the block arm is omitted rather than
#     carried over. Re-add it here if/when window mode is refreshed.
#
#  3. Writes a NEW table (benchmark_table_4tool_chrom.tsv) per the roadmap's
#     "new output dirs, don't overwrite" rule — the old table stays for the delta.
#
# hapFIRE and vg outputs are reused untouched (competitor-side, estimator-independent).
# Usage: bash score_all_competitors_chrom.sh
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate

RES=benchmarks/accuracy_vs_competitors/results
W=benchmarks/accuracy_vs_competitors/work
KM=benchmarks/p80/results/kmate_chrom_p80_filt2inv
TABLE=benchmarks/benchmark_table_4tool_chrom.tsv
SNP_PANEL=$W/shared_snps_p80_Chr1.vcf.gz
META=benchmarks/p80/data/var_pa_p80.meta.npz
CALLED=benchmarks/p80/data/var_pa_p80.var_called.npz
rm -f "$TABLE"

# pool <TAB> refreshed-kMate-tsv. Only pools that have BOTH hapFIRE and vg (SNP+SV)
# outputs are scoreable into the 4-tool table; that is exactly this pair.
POOLS="
cov10_n80_g0_s42_hotspots_p80_chr1|$KM/n80_g0/p80_chrom_filt2inv_n80_g0_cov10_s42.tsv
cov50_n80_g0_s42_hotspots_p80_chr1|$KM/n80_g0/p80_chrom_filt2inv_n80_g0_cov50_s42.tsv
"

done=0; skip=0
for entry in $POOLS; do
  p=${entry%%|*}; KG=${entry##*|}
  HF=$RES/hapfire_${p}_snp_frequency.txt
  VS=$RES/vg_${p}_snp_frequency.txt
  VSV=$W/vg_sv_${p}.tsv
  miss=""
  for f in "$KG" "$HF" "$VS" "$VSV"; do [ -s "$f" ] || miss="$miss $(basename $f)"; done
  if [ -n "$miss" ]; then echo "SKIP $p (missing:$miss)"; skip=$((skip+1)); continue; fi
  echo "== $p"
  python benchmarks/accuracy_vs_competitors/scripts/build_4tool_table.py \
    --pool "$p" --truth benchmarks/p80/sims/$p/recomb_truth.tsv.gz \
    --snp-panel "$SNP_PANEL" --var-meta "$META" --var-called "$CALLED" \
    --kmate-global "$KG" --kmate-mode chrom \
    --hapfire "$HF" --vg-snp "$VS" --vg-sv "$VSV" \
    --table "$TABLE"
  done=$((done+1))
done
echo "scored $done pools, skipped $skip -> $TABLE"
