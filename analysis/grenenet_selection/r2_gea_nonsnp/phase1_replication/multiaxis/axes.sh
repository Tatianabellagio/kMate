# Shared axis/class definitions for the multi-axis clq0.9 GEA extension.
# 20 climate axes (bio1..bio19 + pc1) x 3 variant classes (snp, sv, smallindel) = 60 tasks.
# SLURM array index i (0..59): axis = AXES[i/3], cls = CLS[i%3].
# (Consolidated 2026-07-21: the pooled `nonsnp` class was retired for the 3-class split.)
AXES=(bio1 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 pc1)
CLS=(snp sv smallindel)
axis_for() { echo "${AXES[$(( $1 / 3 ))]}"; }
cls_for()  { echo "${CLS[$(( $1 % 3 ))]}"; }
