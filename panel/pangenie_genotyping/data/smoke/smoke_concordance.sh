#!/bin/bash
# =============================================================================
# smoke_concordance.sh — synthetic-VCF smoke test for loo_concordance.py.
#
# Builds tiny truth + pangenie VCFs (5 records, one per size_class), runs the
# concordance script, and asserts the expected GC / nRD / counts. Designed to
# fail loudly on any silent-output bug before LOO real-data lands.
#
#   chr  pos     class        truth   pg          expected
#   1    1000    SNP          0/0     0/0         match (no nonref)
#   1    2000    small_indel  0/1     0/1         match (nonref both, GC=1)
#   1    3000    small_sv     1/1     0/1         mismatch (nonref both, disc)
#   1    4000    medium_sv    0/0     ./.         pg_missing
#   1    5000    large_sv     0/1     1/1         mismatch (nonref both, disc)
# =============================================================================
set -eo pipefail
# pang env has bgzip + tabix + pysam. Activate before set -u (PYTHONPATH hook).
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pang
set -u

SMOKE=/global/scratch/users/tbellg/kmate/panel/pangenie_genotyping/data/smoke
SCRIPT=/global/scratch/users/tbellg/kmate/panel/pangenie_genotyping/scripts/loo_concordance.py
mkdir -p $SMOKE
cd $SMOKE

# The medium_sv ALT is exactly 1000 chars (max(ref,alt)=1000 → medium_sv class).
# The large_sv ALT is exactly 6000 chars (max=6000 → large_sv class).
MED_ALT=$(python -c 'print("A"*1000)')
LRG_ALT=$(python -c 'print("A"*6000)')

cat > truth.vcf <<EOF
##fileformat=VCFv4.2
##contig=<ID=1,length=10000000>
##INFO=<ID=.,Number=0,Type=Flag,Description="placeholder">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TRUTH
1	1000	rs_snp	A	G	.	PASS	.	GT	0/0
1	2000	rs_ind	A	ACGT	.	PASS	.	GT	0/1
1	3000	rs_ssv	A	$(python -c 'print("C"*100)')	.	PASS	.	GT	1/1
1	4000	rs_msv	A	$MED_ALT	.	PASS	.	GT	0/0
1	5000	rs_lsv	A	$LRG_ALT	.	PASS	.	GT	0/1
EOF

cat > pang.vcf <<EOF
##fileformat=VCFv4.2
##contig=<ID=1,length=10000000>
##INFO=<ID=.,Number=0,Type=Flag,Description="placeholder">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	PG
1	1000	rs_snp	A	G	.	PASS	.	GT	0/0
1	2000	rs_ind	A	ACGT	.	PASS	.	GT	0/1
1	3000	rs_ssv	A	$(python -c 'print("C"*100)')	.	PASS	.	GT	0/1
1	4000	rs_msv	A	$MED_ALT	.	PASS	.	GT	./.
1	5000	rs_lsv	A	$LRG_ALT	.	PASS	.	GT	1/1
EOF

bgzip -f truth.vcf && tabix -fp vcf truth.vcf.gz
bgzip -f pang.vcf  && tabix -fp vcf pang.vcf.gz

python $SCRIPT \
  --pangenie-vcf pang.vcf.gz \
  --truth-vcf    truth.vcf.gz \
  --sample       PG \
  --truth-sample TRUTH \
  --out          smoke_out

echo
echo "=== Asserting expected per-class metrics ==="
# Expected (size_class, n_total, n_matched, n_concordant, GC, n_nonref_either, n_nonref_disc, nRD, n_truth_missing, n_pg_missing)
EXPECT=$(cat <<EOF
SNP	1	1	1	1.0000	0	0	nan	0	0
small_indel	1	1	1	1.0000	1	0	0.0000	0	0
small_sv	1	1	0	0.0000	1	1	1.0000	0	0
medium_sv	1	0	0	nan	0	0	nan	0	1
large_sv	1	1	0	0.0000	1	1	1.0000	0	0
EOF
)

GOT=$(tail -n +2 smoke_out_summary.tsv)
echo "EXPECT:"; echo "$EXPECT"
echo "GOT:";    echo "$GOT"

if [ "$EXPECT" = "$GOT" ]; then
    echo "SMOKE TEST PASSED"
    exit 0
else
    echo "SMOKE TEST FAILED — output != expected"
    diff <(echo "$EXPECT") <(echo "$GOT") || true
    exit 1
fi
