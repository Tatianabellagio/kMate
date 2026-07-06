#!/bin/bash
#SBATCH --job-name=lfmm_hap
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=3:00:00
#SBATCH --output=analysis/grenenet_gea/phase1_replication/logs/lfmm_hap_%j.out
#SBATCH --error=analysis/grenenet_gea/phase1_replication/logs/lfmm_hap_%j.out
# Structure-corrected LFMM (K=16) on the clq0.9 haploblock last-gen matrix (raw hap freq,
# same Y as the Kendall run) -> calibrated per-haploblock p -> re-run the size-matched
# SV-enrichment on the calibrated p. Compares lambda vs the raw Kendall (8.95).
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p analysis/grenenet_gea/phase1_replication/logs results/grenenet_gea/phase1_replication/lfmm
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PR=analysis/grenenet_gea/phase1_replication
LDIR=results/grenenet_gea/phase1_replication/lfmm
CMDIR=results/grenenet_gea/phase1_replication/class_matrices
K=16

echo "== build hap LFMM input (raw freq, z bio1) =="
$PY -u $PR/build_lfmm_input_hap.py --climate bio1

echo "== LFMM ridge K=$K + gif =="
$R $PR/run_lfmm_lastgen.R "$LDIR/lfmm_hap_gen9" $K "$LDIR/lfmm_hap_gen9_bio1_calibp.csv"

echo "== join calibrated p -> per-haploblock GEA csv (hap_id-keyed) =="
$PY - <<'EOF'
import pandas as pd
CMDIR="results/grenenet_gea/phase1_replication/class_matrices"
LDIR="results/grenenet_gea/phase1_replication/lfmm"
rec=pd.read_csv(f"{CMDIR}/hap_gen9.records.csv")            # hap_id-ordered == Y columns
pv=pd.read_csv(f"{LDIR}/lfmm_hap_gen9_bio1_calibp.csv")["pval"].to_numpy()
assert len(pv)==len(rec),(len(pv),len(rec))
rec=rec.assign(pval=pv).rename(columns={"maf":"MAF"})
out=f"{LDIR}/lfmm_hap_gen9_bio1.csv"
rec[["hap_id","chrom","pos","ref_len","alt_len","MAF","block","pval"]].to_csv(out,index=False)
print("  wrote",out,"| min p=%.2e | median p=%.3f"%(rec.pval.min(),rec.pval.median()))
EOF

echo "== size-matched SV-enrichment on LFMM calibrated p =="
$PY -u $PR/hap_sv_enrichment.py \
    --gea "$LDIR/lfmm_hap_gen9_bio1.csv" \
    --out results/grenenet_gea/phase1_replication/hap_sv_enrichment_lfmm.csv
echo "== ALLDONE_LFMM_HAP =="
