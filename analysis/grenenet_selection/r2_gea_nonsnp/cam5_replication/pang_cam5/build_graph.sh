#!/bin/bash
# Build a small founder pangenome graph of CAM5 from the production panel VCF.
# Reference path + variant bubbles (SNP/indel/SV) for the 231 founders.
set -euo pipefail
VG=/global/home/users/tbellg/miniforge3/envs/pangraph/bin
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa
cd "$(dirname "$0")"

# 1) full-gene graph -------------------------------------------------------
$VG/vg construct -r "$REF" -v cam5_region.vcf.gz -R Chr2:11531800-11534533 -a -f > cam5.vg
$VG/vg stats -z -l cam5.vg | tee cam5_graph_stats.txt
$VG/vg view -g cam5.vg > cam5.gfa
# odgi 1D overview (linear pangenome view)
$VG/odgi build -g cam5.gfa -o cam5.og
$VG/odgi viz -i cam5.og -o cam5_odgi_viz.png -x 1600 -y 200

# 2) focused graph of the GEA haplotype (3' end, ~11,533,880-11,534,160) ----
$VG/vg construct -r "$REF" -v cam5_region.vcf.gz -R Chr2:11533880-11534160 -a -f > cam5_hap.vg
# readable bubble graph via graphviz (small region)
$VG/vg view -dpn cam5_hap.vg | $VG/dot -Tpng -o cam5_hap_graph.png
$VG/vg view -dpn cam5_hap.vg | $VG/dot -Tsvg -o cam5_hap_graph.svg
echo "DONE: cam5_odgi_viz.png (overview), cam5_hap_graph.png (3' GEA-haplotype bubbles)"
