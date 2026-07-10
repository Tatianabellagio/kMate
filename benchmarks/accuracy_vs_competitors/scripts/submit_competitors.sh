#!/bin/bash
# Submit the competitor (hapFIRE + vg-SNP + vg-SV) campaign over all p80 benchmark
# pools. Gates on existing outputs so re-runs only fill gaps. p80 ONLY (no p231 vg
# graph; paper decision to keep the competitor comparison on the home-turf panel).
# Usage: bash submit_competitors.sh [--dry]
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
MODE=${1:-both}   # snp | sv | both
DRY=${2:-}
RES=benchmarks/accuracy_vs_competitors/results
W=benchmarks/accuracy_vs_competitors/work
SC=benchmarks/accuracy_vs_competitors/scripts

# canonical p80 grid = pools kMate was scored on (benchmark_table.tsv)
POOLS=$(tail -n +2 benchmarks/benchmark_table.tsv | awk -F'\t' '{print $NF}' \
        | sed -E 's/_(global|block_dynldK500|block)\.tsv$//' | sort -u | grep p80)

sub() { [ "$DRY" = "--dry" ] && { echo "DRY $*"; echo "DRY"; } || sbatch --parsable "$@"; }

n_snp=0; n_sv=0
for p in $POOLS; do
  [ -f benchmarks/p80/sims/$p/reads/r1.fq ] || { echo "skip (no reads): $p"; continue; }
  # SNP-side: need hapFIRE + vg-SNP
  if [ "$MODE" != sv ] && { [ ! -s "$RES/hapfire_${p}_snp_frequency.txt" ] || [ ! -s "$RES/vg_${p}_snp_frequency.txt" ]; }; then
    jid=$(sub $SC/run_compet_snp.sbatch "$p"); echo "SNP  $jid  $p"; n_snp=$((n_snp+1))
  fi
  # SV-side: need vg-SV svidx,est
  if [ "$MODE" != snp ] && [ ! -s "$W/vg_sv_${p}.tsv" ]; then
    jid=$(sub $SC/run_compet_sv.sbatch "$p"); echo "SV   $jid  $p"; n_sv=$((n_sv+1))
  fi
done
echo "submitted: $n_snp SNP-side + $n_sv SV-side jobs"
