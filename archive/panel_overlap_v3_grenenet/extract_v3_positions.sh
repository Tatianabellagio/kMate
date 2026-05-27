#!/bin/bash -l
#SBATCH --job-name=v3_overlap
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/panel_overlap_v3_grenenet/logs/extract_v3_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/panel_overlap_v3_grenenet/logs/extract_v3_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
# =============================================================================
# Extract v3 panel position lists for overlap analysis vs GrENE-Net 231 SNP catalog.
#
# v3 panel = founders_231_chr.haploid.vcf.gz (80 cactus + 151 PG, all 231
# GrENE-Net founders, SNPs+indels+SVs in one file). Same 231 founders as
# GrENE-Net so NO subset complication needed — direct position overlap.
#
# Three position lists:
#   v3_all_positions.tsv         — all (chrom, pos) in v3 (incl. multi-allelic)
#   v3_snp_positions.tsv         — biallelic SNP positions only
#   v3_snp_positions_dedup.tsv   — same, dedup'd to unique (chrom, pos)
# Plus chrom-name-normalized versions (Chr-prefix stripped to match GrENE-Net).
# =============================================================================
set -euo pipefail
eval "$(conda shell.bash hook)" 2>/dev/null || true
BCF=/home/tbellagio/miniforge3/envs/BIOS424/bin/bcftools

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/panel_overlap_v3_grenenet
DATA=$BASE/data
V3_VCF=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/data/merged/founders_231_chr.haploid.vcf.gz

mkdir -p $DATA

# 1) All v3 positions (any variant type, any allele)
echo "[$(date)] dumping ALL v3 positions"
$BCF query -f '%CHROM\t%POS\n' "$V3_VCF" > $DATA/v3_all_positions.tsv
wc -l $DATA/v3_all_positions.tsv

# 2) SNP positions only (biallelic SNPs to keep it apples-to-apples with GrENE-Net)
echo "[$(date)] dumping v3 SNP positions (biallelic only)"
$BCF view -v snps -m2 -M2 "$V3_VCF" -Ou \
  | $BCF query -f '%CHROM\t%POS\n' \
  > $DATA/v3_snp_positions.tsv
wc -l $DATA/v3_snp_positions.tsv

# 3) ALL SNP positions (incl. multi-allelic) — like the old test for ours82
echo "[$(date)] dumping v3 SNP positions (incl. multi-allelic)"
$BCF view -v snps "$V3_VCF" -Ou \
  | $BCF query -f '%CHROM\t%POS\n' \
  | sort -u > $DATA/v3_snp_positions_dedup.tsv
wc -l $DATA/v3_snp_positions_dedup.tsv

# 4) Normalize chrom names (strip 'Chr' prefix) for direct comparison with GrENE-Net
echo "[$(date)] stripping Chr prefix from positions"
for f in v3_all_positions v3_snp_positions v3_snp_positions_dedup; do
  sed 's/^Chr//' $DATA/${f}.tsv > $DATA/${f}.nochr.tsv
done

echo "[$(date)] DONE — file summary:"
ls -la $DATA/v3_*.tsv
