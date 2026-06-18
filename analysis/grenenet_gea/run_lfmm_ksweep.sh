#!/bin/bash
#SBATCH --job-name=lfmm_ksweep
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=8:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/ksweep_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/ksweep_%j.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/analysis/grenenet_gea
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
STEM=/global/scratch/users/tbellg/kmate/results/grenenet_gea/phase1_replication/lfmm/lfmm_snp_gen9
OUT=/global/scratch/users/tbellg/kmate/results/grenenet_gea/lfmm/ksweep_gif_snp_gen9_bio1.csv
mkdir -p "$(dirname "$OUT")"
# SNP-based calibration (apply chosen K to all classes). K=16 = phase-1 value, in the sweep.
$R run_lfmm_ksweep.R "$STEM" "$OUT" "1,2,3,4,6,8,10,12,14,16,18,20"
echo DONE_KSWEEP
