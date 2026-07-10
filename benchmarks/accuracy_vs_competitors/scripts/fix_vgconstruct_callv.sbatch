#!/bin/bash
#SBATCH --job-name=vgcv_fix
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=1:30:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
set -eo pipefail
# Fix the construct-method vg-SV: vg index -x dropped the _alt_ paths (no -L), so
# `vg call -v` found 0. Rebuild the xg WITH -L from the EXISTING .vg (nodes unchanged
# => existing pack still valid, no remap), then call -v the panel SV VCF.
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vg
T=${SLURM_CPUS_PER_TASK:-16}
POOL=cov10_n80_g0_s42_hotspots_p80_chr1
W=benchmarks/accuracy_vs_competitors/work
VG=$W/c_${POOL}.vg
PACK=$W/${POOL}_construct.pack
SVVCF=$W/panel_sv_p80_Chr1.vcf.gz
XGL=$W/c_${POOL}.alts.xg
OUT=$W/${POOL}_vgconstruct_sv_FIXED.vcf

echo "[$(date)] host=$(hostname) rebuild xg WITH -L (keep _alt_ paths)"
vg index -x "$XGL" -L "$VG" -t "$T"
echo "[$(date)] xg done: alt paths now = $(vg paths -L -x "$XGL" | grep -c '_alt_')"
echo "[$(date)] vg call -v (panel SVs) on existing pack"
vg call "$XGL" -k "$PACK" -v "$SVVCF" -t "$T" > "$OUT"
echo "[$(date)] DONE: $(grep -vc '^#' "$OUT") variant lines -> $OUT"
echo "=== sample records w/ AD ==="
grep -v '^#' "$OUT" | awk -F'\t' 'length($5)-length($4)>=50 || length($4)-length($5)>=50' | head -3 | cut -f1-5,10
