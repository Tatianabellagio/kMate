#!/bin/bash
#SBATCH --job-name=arch3_annot
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/annotate_%j.out
#SBATCH --error=logs/annotate_%j.err

mkdir -p logs
set -euo pipefail

cd /global/scratch/users/tbellg/kmate/scratch/arch3_test

GFA=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output/pang_1001gplus_all.gfa.gz
SCRIPT=/global/scratch/users/tbellg/kmate/external_tools/genotyping-pipelines/prepare-vcf-MC/workflow/scripts/annotate_vcf.py
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

VCF=cactus_78_test_chr1_5_14M.vcf
OUTPREFIX=cactus_78_test_annotated

echo "[$(date)] Starting annotate_vcf on test region Chr1:5.8M-14M..."
echo "VCF: $VCF ($(ls -la $VCF | awk '{print $5}') bytes)"
echo "GFA: $GFA ($(ls -la $GFA | awk '{print $5}') bytes)"

$PY -u $SCRIPT -vcf $VCF -gfa $GFA -o $OUTPREFIX

echo "[$(date)] annotate_vcf complete."
echo
echo "=== Output files ==="
ls -la ${OUTPREFIX}*.vcf

echo
echo "=== Multi-allelic output: header check ==="
grep "^##INFO=<ID=ID" ${OUTPREFIX}.vcf | head -2

echo
echo "=== Multi-allelic: 3 spot-check bubbles ==="
for pos in 5869846 10421232 13843898; do
  echo "--- Chr1:$pos ---"
  grep -P "^Chr1\t$pos\t" ${OUTPREFIX}.vcf | head -1 | awk -F'\t' '{
    print "  REF_len="length($4)" N_alts="gsub(",","&",$5)+1;
    n=split($8,a,";");
    for(i=1;i<=n;i++) if(a[i] ~ /^ID=/) {print "  INFO/ID first 200 chars: " substr(a[i],4,200); break}
  }'
done

echo
echo "=== Biallelic catalog: count and spot-check ==="
echo "Total biallelic records: $(grep -v '^#' ${OUTPREFIX}_biallelic.vcf | wc -l)"
echo
echo "Biallelic records at exactly pos 5870018, 10421645, 13843898:"
for pos in 5870018 10421645 13843898; do
  n=$(awk -v p=$pos 'BEGIN{n=0} !/^#/ && $2==p {n++} END{print n}' ${OUTPREFIX}_biallelic.vcf)
  echo "  pos $pos: $n biallelic records"
  awk -v p=$pos '!/^#/ && $2==p {print "    REF="$4" ALT="$5" ID-col3="$3" INFO="$8}' ${OUTPREFIX}_biallelic.vcf | head -3
done
