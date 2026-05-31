#!/bin/bash
#SBATCH --job-name=arch3_build
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --array=2-5
#SBATCH --output=logs/arch3_build_%A_%a.out
#SBATCH --error=logs/arch3_build_%A_%a.err
#
# Generic arch3 A1->A5 driver, one array task per chromosome.
# arch3 is THE sole production decomposition (see docs/PIPELINE_STATE.md §0).
# Materializes per-chrom job scripts from the validated chr1/ templates by
# substituting chr1->chrN, Chr1->ChrN, and the removed conda envs -> kmate,
# then runs A1 (annotate) -> A2 (PG) -> A3 (cactus) -> A4 (merge) -> A5 (var_pa).
#
# Output per chrom: panel/arch3/chr{N}/merged_231_chr{N}_final.vcf.gz
#                   panel/arch3/chr{N}/var_pa_231_arch3_chr{N}.{var_pa,var_called,meta}.npz
#
# Default array is 2-5 (Chr1 already built + validated). Override:
#   sbatch --array=2 panel/arch3/build_arch3_chrom_array.sh        # one chrom
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/kmate
ARCH=$BASE/panel/arch3
N=${SLURM_ARRAY_TASK_ID:-2}
CHR=Chr${N}; chr=chr${N}
SRC=$ARCH/chr1
DST=$ARCH/${chr}
mkdir -p "$DST/logs"

echo "[$(date)] === arch3 build $CHR -> $DST ==="

# 1) Materialize per-chrom job scripts from the chr1 templates (A1..A5, globbed
#    so we don't assume the _chr1 suffix — A5's template is jobA5_build_cnvar.sh).
#    s/chr1/chrN/  -> filenames + output dir   (lowercase var values)
#    s/Chr1/ChrN/  -> REGION / bcftools -r
#    env remap     -> kmate (gwas/hapfm/sequencing_pipeline/pang/pangenie/basic removed on this cluster)
#    NOTE: uppercase token "CHR1" in internal var NAMES (VCF_CHR1, PG_CHR1, ...) is
#    intentionally left as-is — it's a cosmetic variable name, internally consistent.
GEN=()
for srcf in $(ls $SRC/jobA[1-5]_*.sh | sort); do
    dstf=$DST/$(basename "$srcf" | sed "s/chr1/${chr}/")
    sed -e "s/chr1/${chr}/g" -e "s/Chr1/${CHR}/g" \
        -e "s#miniforge3/envs/\(gwas\|hapfm\|sequencing_pipeline\|pang\|pangenie\|basic\)/bin#miniforge3/envs/kmate/bin#g" \
        "$srcf" > "$dstf"
    chmod +x "$dstf"
    GEN+=("$dstf")
done
[ ${#GEN[@]} -eq 5 ] || { echo "ERROR: expected 5 stage scripts, got ${#GEN[@]}"; exit 1; }

# 2) Guard: no stale lowercase chr1 / mixed Chr1 / dead-env token should remain.
if grep -nE "chr1|Chr1|envs/(gwas|hapfm|sequencing_pipeline|pang|pangenie|basic)/bin" "${GEN[@]}"; then
    echo "ERROR: stale chr1/Chr1/dead-env token left after templating (see above)"; exit 1
fi
echo "[$(date)] templating clean: ${#GEN[@]} stage scripts in $DST"

# 3) Run the stages in dependency order (plain bash; each cd's into $DST itself).
for f in "${GEN[@]}"; do
    echo; echo "[$(date)] ===== RUN $(basename "$f") ($CHR) ====="
    bash "$f"
done

echo
echo "[$(date)] === arch3 $CHR COMPLETE ==="
ls -lh $DST/merged_231_${chr}_final.vcf.gz $DST/var_pa_231_arch3_${chr}.*.npz
