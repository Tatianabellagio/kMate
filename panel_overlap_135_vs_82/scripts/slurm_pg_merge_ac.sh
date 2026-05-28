#!/bin/bash
#SBATCH --job-name=overlap_pg_merge
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/pg_merge_%j.out
#SBATCH --error=logs/pg_merge_%j.err
mkdir -p logs
set -euo pipefail

# Merge the 151 PG-genotyped per-founder VCFs (from pangenie_genotyping/data/genotyped/),
# fill AC/AN over the 151, then emit per-ALT AC at the same alt-idx granularity
# as Tier 1 ac_82 / ac_53 tables.
#
# Output: panel_overlap_135_vs_82/results/ac_pg151_Chr{1..5}.tsv.gz
#  (schema matches ac_82 / ac_53 tables: chrom pos ref_len n_alt alt_idx alt_len ac an)
#
# Notes:
# - 151 sources: pangenie_genotyping/data/genotyped/*_genotyping.vcf.gz
#   (matches missing_151_ecotypes.txt exactly; these are the GrENE-Net founders
#   that lack long-read assemblies and were genotyped via PanGenie against pang_135).
# - We do NOT apply v3qc het-mask or V4 filter here: this is the raw PG signal,
#   which is what the question asks about ("how much information does PG provide
#   at extras-introduced bubbles before any QC").

source /global/home/users/tbellg/miniforge3/etc/profile.d/conda.sh
conda activate sequencing_pipeline

BASE=/global/scratch/users/tbellg/hapfire_sv
WORK=$BASE/panel_overlap_135_vs_82
GENOTYPED=$BASE/pangenie_genotyping/data/genotyped

echo "[$(date +%H:%M:%S)] discovering input PG VCFs"

# 151 from genotyped/
ls $GENOTYPED/*_genotyping.vcf.gz 2>/dev/null | sort > $WORK/data/pg_vcfs.txt
N_TOTAL=$(wc -l < $WORK/data/pg_vcfs.txt)
echo "  total PG VCFs to merge: $N_TOTAL"

# Sanity check: all VCFs must be indexed
echo "[$(date +%H:%M:%S)] checking indices"
while read f; do
  [ -s "${f}.tbi" ] || tabix -p vcf -f "$f"
done < $WORK/data/pg_vcfs.txt

# Merge per chrom (cheaper than one whole-genome merge)
mkdir -p $WORK/results
for CHROM in Chr1 Chr2 Chr3 Chr4 Chr5; do
  out=$WORK/results/ac_pg151_${CHROM}.tsv.gz
  if [ -s "$out" ]; then
    echo "[$(date +%H:%M:%S)] $CHROM: $out exists, skipping"
    continue
  fi
  echo "[$(date +%H:%M:%S)] merging $CHROM"
  bcftools merge --threads 6 -l $WORK/data/pg_vcfs.txt -r $CHROM -m none -Ou \
    | bcftools +fill-tags --threads 2 -- -t AC,AN \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' \
    | awk -v OFS='\t' '
        BEGIN {print "chrom\tpos\tref_len\tn_alt\talt_idx\talt_len\tac_pg\tan_pg"}
        {
          n_alt = split($4, alts, ",");
          split($5, acs, ",");
          ref_len = length($3);
          an = $6;
          for (i=1; i<=n_alt; i++) {
            alt_len = length(alts[i]);
            ac = (i in acs) ? acs[i] : 0;
            print $1, $2, ref_len, n_alt, i, alt_len, ac, an;
          }
        }' \
    | gzip > $out
  echo "[$(date +%H:%M:%S)] $CHROM: $(zcat $out | wc -l) rows → $out"
done

echo "[$(date +%H:%M:%S)] all PG chroms done"
