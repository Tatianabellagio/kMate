#!/bin/bash
#SBATCH --job-name=nonsnp_ksweep
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=8:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/nonsnp_ksweep_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/nonsnp_ksweep_%j.out
# STAGE A: build the pooled non-SNP (indel+SV) class matrix + LFMM input, then
# run the GIF K-sweep on it (mirrors the existing SNP sweep) so we can pick a
# common, well-calibrated K for the SNP-vs-nonSNP LFMM comparison.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PR=analysis/grenenet_gea/phase1_replication
LDIR=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/results/lfmm
KOUT=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/lfmm/ksweep_gif_nonsnp_gen9_bio1.csv
mkdir -p "$LDIR" "$(dirname "$KOUT")"

echo "== build pooled non-SNP class matrix (gen9, MAF>=0.05) =="
$PY -u $PR/build_class_matrices.py --classes nonsnp --gens 9

echo "== build LFMM input (Δp, impute, env) for nonsnp gen9 =="
$PY -u $PR/build_lfmm_input.py --class nonsnp --gen 9 --climate bio1

echo "== GIF K-sweep on nonsnp =="
$R analysis/grenenet_gea/run_lfmm_ksweep.R \
   "$LDIR/lfmm_nonsnp_gen9" "$KOUT" "1,2,4,6,8,10,12,14,16,18,20"
echo DONE_NONSNP_KSWEEP
