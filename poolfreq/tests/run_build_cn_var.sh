#!/bin/bash
#SBATCH --job-name=build_cn_var
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/build_cn_var_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/build_cn_var_%j.err

set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq

# Use the biallelic.norm VCF — same one used as PanGenie input
# (one row per biallelic ALT of each top-level bubble)
/home/tbellagio/miniforge3/envs/hapfm/bin/python src/build_cn_var.py \
    --vcf /carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_test/raw_vcfbub_lv0_diploid.biallelic.norm.vcf.gz \
    --out data/cn_var_82
