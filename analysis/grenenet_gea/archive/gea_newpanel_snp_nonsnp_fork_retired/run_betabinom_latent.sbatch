#!/bin/bash
#SBATCH --job-name=betabinom_lf
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=analysis/grenenet_gea/gea_newpanel/betabinom_lf_%j.out
#SBATCH --error=analysis/grenenet_gea/gea_newpanel/betabinom_lf_%j.out
# Beta-binomial (overdispersion-aware) + K=16 latent-factor climate-GEA, gen9, bio1.
# Fixed-rho beta-binomial == effective-sample-size binomial (see script docstring).
# Stage 1: Williams moment rho (K=16 residuals, MAF>=0.05 subsample) per class.
# Stage 2: genome pass BB K=16 (primary) and BB K=0 (structure-vs-overdisp decomp).
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
G=analysis/grenenet_gea/gea_newpanel
T="${SLURM_CPUS_PER_TASK:-16}"
echo "== host $(hostname) | beta-binomial K=16 GEA gen9 bio1 =="

for CLS in nonsnp snp; do
    echo "---- rho estimation: $CLS ----"
    $PY -u $G/run_betabinom_latent.py --class "$CLS" --estimate-rho \
        --nsub 15000 --rho-maf 0.05 --threads "$T"
done

for CLS in nonsnp snp; do
    for K in 16 0; do
        echo "---- genome pass: $CLS K=$K ----"
        $PY -u $G/run_betabinom_latent.py --class "$CLS" --gen 9 --climate bio1 \
            --k "$K" --rho-json "rho_${CLS}.json" --threads "$T"
    done
done
echo "== ALLDONE_BETABINOM_LF =="
