#!/bin/bash
# Tier 1 (fast): per-ALT AC over the 82 GrENE-Net cactus columns vs the 53
# extras columns in the pang_135 raw VCF.
#
# Key trick: multi-allelics are NOT decomposed. AC has Number=A (per-ALT) and
# PanGenie inherits the panel's ALT ordering. The join key with PG VCFs is
# (chrom, pos, alt_idx), not the literal ALT string — which saves us from
# emitting/comparing 100kb+ SV strings.
#
# Output schema (one row per (record, ALT_i)):
#   chrom  pos  ref_len  n_alt  alt_idx  alt_len  ac_82  an_82  ac_53  an_53
set -euo pipefail

source /global/home/users/tbellg/miniforge3/etc/profile.d/conda.sh
conda activate sequencing_pipeline

WORK=/global/scratch/users/tbellg/kmate/panel_overlap_135_vs_82
RAW=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.raw.vcf.gz
GROUP82=$WORK/data/grenenet_82_in_pang135.txt
GROUP53=$WORK/data/extras_53.txt
OUTDIR=$WORK/results
CHROM=${1:-Chr1}

echo "[$(date +%H:%M:%S)] Tier 1 AC split for $CHROM (fast path, no norm)"

make_ac_table() {
  local group_file=$1
  local tag=$2
  local out=$3
  echo "[$(date +%H:%M:%S)]   ac_${tag} for $CHROM ..."
  bcftools view --threads 2 -r "$CHROM" -S "$group_file" --force-samples "$RAW" \
    | bcftools +fill-tags --threads 2 -- -t AC,AN \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' \
    | awk -v OFS='\t' -v tag="$tag" '
        BEGIN {print "chrom\tpos\tref_len\tn_alt\talt_idx\talt_len\tac_"tag"\tan_"tag}
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
    | gzip > "$out"
  echo "[$(date +%H:%M:%S)]   wrote $out ($(zcat $out | wc -l) lines)"
}

make_ac_table "$GROUP82" 82 "$OUTDIR/ac_82_${CHROM}.tsv.gz" &
PID82=$!
make_ac_table "$GROUP53" 53 "$OUTDIR/ac_53_${CHROM}.tsv.gz" &
PID53=$!

wait $PID82 $PID53

echo "[$(date +%H:%M:%S)] joining 82 + 53 for $CHROM"

paste <(zcat "$OUTDIR/ac_82_${CHROM}.tsv.gz") \
      <(zcat "$OUTDIR/ac_53_${CHROM}.tsv.gz") \
  | awk -v OFS='\t' '
      NR==1 {print "chrom\tpos\tref_len\tn_alt\talt_idx\talt_len\tac_82\tan_82\tac_53\tan_53"; next}
      {
        if ($1!=$9 || $2!=$10 || $5!=$13) {
          print "MISMATCH at line "NR": "$1":"$2"."$5" vs "$9":"$10"."$13 > "/dev/stderr"; exit 2
        }
        print $1,$2,$3,$4,$5,$6,$7,$8,$15,$16
      }' \
  | gzip > "$OUTDIR/ac_split_${CHROM}.tsv.gz"

echo "[$(date +%H:%M:%S)] done. Output: $OUTDIR/ac_split_${CHROM}.tsv.gz"
echo "[$(date +%H:%M:%S)] rows: $(zcat $OUTDIR/ac_split_${CHROM}.tsv.gz | wc -l)"
echo "[$(date +%H:%M:%S)] head:"
zcat "$OUTDIR/ac_split_${CHROM}.tsv.gz" | head -5
