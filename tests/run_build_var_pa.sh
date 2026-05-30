#!/bin/bash
#SBATCH --job-name=build_var_pa
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/build_var_pa_%j.out
#SBATCH --error=logs/build_var_pa_%j.err

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate

# Use the biallelic.norm VCF — same one used as PanGenie input
# (one row per biallelic ALT of each top-level bubble)
/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python src/build_var_pa.py \
    --vcf /global/scratch/users/tbellg/kmate/pangenie_test/raw_vcfbub_lv0_diploid.biallelic.norm.vcf.gz \
    --out data/var_pa_82
