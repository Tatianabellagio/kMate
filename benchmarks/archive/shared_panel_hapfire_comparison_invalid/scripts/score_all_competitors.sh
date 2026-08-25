#!/bin/bash
# Score every p80 pool through build_4tool_table.py into the combined 4-tool table.
# Gates on all required estimates being present (kMate g/b + hapFIRE + vg-SNP + vg-SV).
# Rebuilds the table from scratch each run (idempotent). Usage: bash score_all_competitors.sh
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
RES=benchmarks/accuracy_vs_competitors/results
W=benchmarks/accuracy_vs_competitors/work
TSV=benchmarks/benchmark_runs/tsv
TABLE=benchmarks/benchmark_table_4tool.tsv
SNP_PANEL=$W/shared_snps_p80_Chr1.vcf.gz
META=benchmarks/p80/data/var_pa_p80.meta.npz
CALLED=benchmarks/p80/data/var_pa_p80.var_called.npz
rm -f "$TABLE"

POOLS=$(tail -n +2 benchmarks/benchmark_table.tsv | awk -F'\t' '{print $NF}' \
        | sed -E 's/_(global|block_dynldK500|block)\.tsv$//' | sort -u | grep p80)

done=0; skip=0
for p in $POOLS; do
  KG=$TSV/${p}_global.tsv; KB=$TSV/${p}_block_dynldK500.tsv
  HF=$RES/hapfire_${p}_snp_frequency.txt; VS=$RES/vg_${p}_snp_frequency.txt; VSV=$W/vg_sv_${p}.tsv
  miss=""
  for f in "$KG" "$KB" "$HF" "$VS" "$VSV"; do [ -s "$f" ] || miss="$miss $(basename $f)"; done
  if [ -n "$miss" ]; then echo "SKIP $p (missing:$miss)"; skip=$((skip+1)); continue; fi
  python benchmarks/accuracy_vs_competitors/scripts/build_4tool_table.py \
    --pool "$p" --truth benchmarks/p80/sims/$p/recomb_truth.tsv.gz \
    --snp-panel "$SNP_PANEL" --var-meta "$META" --var-called "$CALLED" \
    --kmate-global "$KG" --kmate-block "$KB" --hapfire "$HF" --vg-snp "$VS" --vg-sv "$VSV" \
    --table "$TABLE" 2>&1 | grep -E 'SNP allrec|SV +full' | sed "s/^/  [$p] /" || true
  done=$((done+1))
done
echo "scored $done pools, skipped $skip -> $TABLE"
