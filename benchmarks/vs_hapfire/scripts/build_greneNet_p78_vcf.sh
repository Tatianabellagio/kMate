#!/bin/bash
# hapFIRE's panel for the long-read-founder pools: the GrENE-Net Chr1 VCF cut to the 78 long-read founders
# (benchmarks/p80/data/samples_78_acc.txt), keeping sites still polymorphic among them. The commands recorded in the
# header of work/greneNet_chr1_p78.vcf.gz (built 2026-10-11).  Usage: bash build_greneNet_p78_vcf.sh
set -euo pipefail
BCF=/global/home/users/tbellg/miniforge3/bin/bcftools
B=/global/scratch/projects/fc_moilab/tbellg/kmate/benchmarks
W=$B/vs_hapfire/work
$BCF view -s "$(paste -sd, $B/p80/data/samples_78_acc.txt)" -Ou $W/greneNet_chr1.vcf.gz \
  | $BCF plugin fill-tags -Ou -- -t AC,AN \
  | $BCF view -i 'AC>0 && AC<AN' -Oz -o $W/greneNet_chr1_p78.vcf.gz
$BCF index -t $W/greneNet_chr1_p78.vcf.gz
