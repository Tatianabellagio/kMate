#!/bin/bash
#SBATCH --job-name=bp_temporal_s4
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --output=analysis/grenenet_gea/logs/bp_temporal_%j.out
#SBATCH --error=analysis/grenenet_gea/logs/bp_temporal_%j.err
# BayPass TEMPORAL (founding->evolved) on haploblocks, site 4 smoke test:
#   estimate Omega via the core model (2 chains, on a 3k-locus subset) -> then OUR fast
#   analytic Omega-conditioned C2 contrast on ALL loci. See BAYPASS_TEMPORAL_PLAN.md.
set -euo pipefail
hostname
cd /global/scratch/users/tbellg/kmate
D=results/grenenet_gea/hapfreq/baypass_temporal_site4

# --- Omega estimation (BayPass core model, baypass env) ---
source ~/miniforge3/etc/profile.d/conda.sh
export LD_LIBRARY_PATH=/global/home/users/tbellg/miniforge3/envs/baypass/lib:$LD_LIBRARY_PATH
conda activate /global/home/users/tbellg/miniforge3/envs/baypass

# 3k-locus subset for Omega (from the current 12-pop geno.txt)
python -c "import numpy as np; G=np.loadtxt('$D/geno.txt',dtype=int); r=np.random.RandomState(0); np.savetxt('$D/geno_sub.txt', G[r.choice(len(G),3000,replace=False)], fmt='%d')"

cd $D
for ch in 1 2; do
  echo "=== core chain $ch  $(date +%T) ==="
  ~/bin/g_baypass -npop 12 -gfile geno_sub.txt -outprefix core_c$ch \
      -nthreads 8 -npilot 20 -pilotlength 1000 -nval 2500 -seed $((ch*7919))
done
cp core_c1_mat_omega.out core_mat_omega.out          # chain-1 Omega = the one we use

# --- analytic Omega-conditioned C2 on ALL loci (kmate env) ---
cd /global/scratch/users/tbellg/kmate
conda deactivate
conda activate kmate
SITE=4 python analysis/grenenet_gea/omega_contrast.py

# FMD-style Omega reproducibility between the 2 chains
python -c "
import numpy as np
A=np.loadtxt('$D/core_c1_mat_omega.out'); B=np.loadtxt('$D/core_c2_mat_omega.out')
print('Omega chain1-vs-chain2 Pearson r =', round(float(np.corrcoef(A.ravel(),B.ravel())[0,1]),4))
"
echo "DONE"
