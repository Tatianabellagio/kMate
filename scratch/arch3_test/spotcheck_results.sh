#!/bin/bash
# Spot-check Arch 3 (HPRC prepare-vcf-MC + convert-to-biallelic) output.
# Computes carrier counts at the 3 canonical positions and compares to known truth.
set -euo pipefail

cd /global/scratch/users/tbellg/kmate/scratch/arch3_test

OUT_VCF=cactus_78_per_sample_biallelic.vcf
[ -s "$OUT_VCF" ] || { echo "ERROR: missing $OUT_VCF"; exit 1; }

echo "=== Per-sample biallelic VCF summary ==="
N_HEADER=$(grep -c "^#" $OUT_VCF)
N_REC=$(grep -vc "^#" $OUT_VCF)
N_SAMP=$(grep "^#CHROM" $OUT_VCF | awk '{print NF-9}')
echo "Header lines: $N_HEADER"
echo "Records: $N_REC"
echo "Samples: $N_SAMP (expected 78)"

echo
echo "=== Spot-check carrier counts ==="
echo
for entry in \
  "5870018:T:A:Chr1:5870018 over-counting case — expected AC≈1-2 (xwu truth)" \
  "10421645:T:C:Chr1:10421645 under-counting case — expected AC≈70-75 (cactus side of 130/135 xwu truth)" \
  "13843898:C:T:Chr1:13843898 C→T atomic (TT path carriers only)" \
  "13843898:CT:TC:Chr1:13843898 CT→TC COMPLEX (TC path carriers)" \
  "13843898:CT:TG:Chr1:13843898 CT→TG COMPLEX (TG path carriers)"
do
  IFS=':' read -r pos ref alt label <<< "$entry"
  echo "--- $label ---"
  awk -v p=$pos -v r=$ref -v a=$alt '
    !/^#/ && $2==p && $4==r && $5==a {
      n_rec++;
      n_0=0; n_1=0; n_dot=0;
      for (i=10; i<=NF; i++) {
        split($i, fmt, ":"); gt=fmt[1];
        gsub(/[|\/]/, " ", gt);
        n_a=split(gt, alleles, " ");
        for (j=1; j<=n_a; j++) {
          if (alleles[j]==".") n_dot++;
          else if (alleles[j]=="0") n_0++;
          else n_1++;
        }
      }
      print "  REF="r" ALT="a"  records="n_rec"  ALT alleles="n_1"  REF alleles="n_0"  missing="n_dot;
      # carrier count (samples with ≥1 ALT allele) — count unique samples
      n_carriers=0; n_called=0; n_total_samp=0;
      for (i=10; i<=NF; i++) {
        n_total_samp++;
        split($i, fmt, ":"); gt=fmt[1];
        gsub(/[|\/]/, " ", gt);
        n_a=split(gt, alleles, " ");
        any_alt=0; any_called=0;
        for (j=1; j<=n_a; j++) {
          if (alleles[j] != ".") any_called=1;
          if (alleles[j] != "." && alleles[j] != "0") any_alt=1;
        }
        if (any_alt) n_carriers++;
        if (any_called) n_called++;
      }
      print "  samples="n_total_samp"  carriers="n_carriers"  called="n_called"  AC/AN="n_carriers"/"n_called;
    }
  ' $OUT_VCF
  echo
done

echo "=== OR-merged carriers at coord 13843898 (across all 3 records, any ALT with T at coord) ==="
echo "Note: requires coord-aware aggregation — see if total matches xwu truth ~109/135"
python3 << 'EOF'
import sys
target_pos = 13843898
samples = None
carriers = None
called = None
records_seen = []

with open("cactus_78_per_sample_biallelic.vcf") as f:
    for line in f:
        if line.startswith("##"): continue
        if line.startswith("#CHROM"):
            samples = line.rstrip().split("\t")[9:]
            carriers = [False] * len(samples)
            called = [False] * len(samples)
            continue
        parts = line.rstrip().split("\t")
        if int(parts[1]) != target_pos: continue
        ref = parts[3]; alt = parts[4]
        # Does ALT[0] == 'T' (i.e., this ALT encodes T at coord 13843898)?
        if len(alt) == 0 or alt[0] != "T": continue
        records_seen.append((parts[3], parts[4]))
        fmt = parts[8].split(":")
        gt_idx = fmt.index("GT") if "GT" in fmt else 0
        for i, samp in enumerate(parts[9:]):
            gt = samp.split(":")[gt_idx]
            if gt == "." or gt == ".|." or gt == "./.": continue
            alleles = gt.replace("|", "/").split("/")
            any_alt = any(a not in ("0", ".") for a in alleles)
            any_called = any(a != "." for a in alleles)
            if any_alt: carriers[i] = True
            if any_called: called[i] = True

print(f"  Records contributing (ALT[0]=='T'): {records_seen}")
print(f"  Total samples: {len(samples)}")
print(f"  OR-merged carriers (have T at coord 13843898): {sum(carriers)}")
print(f"  OR-merged called: {sum(called)}")
EOF
