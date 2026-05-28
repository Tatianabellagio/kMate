#!/bin/bash
#SBATCH --job-name=v3qc_diff
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/v3qc_diff_%j.out
#SBATCH --error=logs/v3qc_diff_%j.err

# 1. Build cn_var_231_v3qc.ref_alt.tsv.gz (matches cn_var_v3qc record order)
# 2. Compute v3qc-vs-GN per-record AF diff → scratch/v3qc_vs_gn_per_record_af.npz
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/kmate
BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

VCF=$BASE/pangenie_genotyping/data/v3qc/founders_231_v3qc.haploid.vcf.gz
REFALT=$BASE/data/cn_var_231_v3qc.ref_alt.tsv.gz

if [ ! -s "$REFALT" ]; then
    echo "[$(date)] Step 1: build $REFALT"
    $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\n' $VCF | gzip > $REFALT
    echo "  lines: $(zcat $REFALT | wc -l)"
fi

echo "[$(date)] Step 2: compute v3qc-vs-GN per-record AF diff"
$PY -u $BASE/scratch/compute_per_record_af_v3qc.py

echo "[$(date)] DONE"
ls -lh $BASE/scratch/v3qc_vs_gn_per_record_af.npz
