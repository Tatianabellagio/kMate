#!/bin/bash
# Unified re-score for one pool: rebuild SNP truth (with physical column) and score
# all available tools (kMate, hapFIRE, vg) vs truth in 3 views:
#   _phys      physical truth (missing->REF), ALL sites      [primary]
#   _fullcalled physical truth, n_called==F only             [convention-free]
#   _mar       MAR truth (kMate estimand)                    [for contrast]
# Usage: bash score_pool.sh <POOL> [PANEL]
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
POOL=${1:?pool}; PANEL=${2:-p80}
SC=benchmarks/accuracy_vs_competitors/scripts
W=benchmarks/accuracy_vs_competitors/work
R=benchmarks/accuracy_vs_competitors/results
COV=$(echo "$POOL" | sed -E 's/^cov([0-9]+)_.*/\1/')
REG=$(echo "$POOL" | sed -E 's/^cov[0-9]+_//; s/_s42_.*$//')
F=$(source ~/miniforge3/etc/profile.d/conda.sh; conda activate kmate; bcftools query -l benchmarks/$PANEL/data/pangenome_${PANEL}_chr1.vcf.gz | wc -l)
KM=$(ls benchmarks/$PANEL/results/cactus_em_global_filt2_mb/$REG/*_cov${COV}_*.tsv 2>/dev/null | head -1)
HF=$R/hapfire_${POOL}_snp_frequency.txt
VG=$R/vg_${POOL}_snp_frequency.txt
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
echo "[score_pool] $POOL  cov=$COV reg=$REG F=$F"
python $SC/build_truth.py \
  --panel-vcf benchmarks/$PANEL/data/pangenome_${PANEL}_chr1.vcf.gz \
  --sites-vcf $W/shared_snps_${PANEL}_Chr1.vcf.gz \
  --weights benchmarks/$PANEL/sims/$POOL/pool_weights.tsv \
  --out $W/truth_${POOL}.tsv 2>/dev/null
conda activate basic
args=(--truth $W/truth_${POOL}.tsv --F $F)
[ -s "$KM" ] && args+=(--kmate "$KM")
[ -s "$HF" ] && args+=(--hapfire "$HF")
[ -s "$VG" ] && args+=(--vg "$VG")
echo "  tools: kMate=$([ -s "$KM" ]&&echo y) hapFIRE=$([ -s "$HF" ]&&echo y) vg=$([ -s "$VG" ]&&echo y)"
echo "### PHYSICAL (all sites) ###"
python $SC/score_competitors.py "${args[@]}" --truth-col truth_af_phys --out-prefix $R/${POOL}_phys 2>/dev/null | cat
echo "### FULLY-CALLED (convention-free) ###"
python $SC/score_competitors.py "${args[@]}" --truth-col truth_af_phys --miss-thr 0.0 --out-prefix $R/${POOL}_fullcalled 2>/dev/null | cat
echo "### MAR (kMate estimand) ###"
python $SC/score_competitors.py "${args[@]}" --truth-col truth_af --out-prefix $R/${POOL}_mar 2>/dev/null | cat
