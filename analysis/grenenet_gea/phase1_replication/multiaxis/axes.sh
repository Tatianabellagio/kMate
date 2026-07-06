# Shared axis/class definitions for the multi-axis clq0.9 GEA extension.
# 20 climate axes (bio1..bio19 + pc1) x 2 variant classes (snp, nonsnp) = 40 tasks.
# SLURM array index i (0..39): axis = AXES[i/2], cls = CLS[i%2].
AXES=(bio1 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 pc1)
CLS=(snp nonsnp)
axis_for() { echo "${AXES[$(( $1 / 2 ))]}"; }
cls_for()  { echo "${CLS[$(( $1 % 2 ))]}"; }
