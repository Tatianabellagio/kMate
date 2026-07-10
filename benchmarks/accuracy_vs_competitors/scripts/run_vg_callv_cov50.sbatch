#!/bin/bash
#SBATCH --job-name=vgcv_cov50
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
set -eo pipefail
# Coverage control for the vg-SV comparison: reuse the EXISTING construct graph
# (depends only on REF+panel VCF, not reads) and run the Ash recipe at cov50 so
# we can attribute vg's cov10 underperformance to coverage vs representation.
# map (cov50 reads) -> pack -Q5 -> call -v panel SVs. Score locally afterwards.
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vg
T=${SLURM_CPUS_PER_TASK:-16}

W=benchmarks/accuracy_vs_competitors/work
GRAPH=$W/c_cov10_n80_g0_s42_hotspots_p80_chr1     # reused graph prefix (.xg/.gcsa)
XGL=$W/c_cov10_n80_g0_s42_hotspots_p80_chr1.alts.xg  # -L xg WITH _alt_ paths (for pack/call)
SVVCF=$W/panel_sv_p80_Chr1.vcf.gz
READS=benchmarks/p80/sims/cov50_n80_g0_s42_hotspots_p80_chr1/reads
GAM=$W/cov50_n80_g0_construct.gam
PACK=$W/cov50_n80_g0_construct.pack
OUT=$W/cov50_n80_g0_s42_hotspots_p80_chr1_vgconstruct_sv_FIXED.vcf

# Stage the big mmap'd indexes to node-local disk (Lustre-mmap stall lesson).
TMP=${TMPDIR:-/tmp}/vgcov50.$SLURM_JOB_ID
mkdir -p "$TMP"
echo "[$(date)] host=$(hostname) staging indexes -> $TMP"
cp "$GRAPH.xg" "$GRAPH.gcsa" "$GRAPH.gcsa.lcp" "$XGL" "$TMP"/ 2>/dev/null || \
  cp "$GRAPH.xg" "$GRAPH.gcsa" "$XGL" "$TMP"/   # .lcp may or may not exist
LGRAPH=$TMP/$(basename "$GRAPH")
LXGL=$TMP/$(basename "$XGL")

echo "[$(date)] STEP map (cov50)"
vg map -d "$LGRAPH" -f "$READS/r1.fq" -f "$READS/r2.fq" -t "$T" > "$GAM"
echo "[$(date)] map done: $(ls -la $GAM | awk '{print $5}') bytes"

echo "[$(date)] STEP pack -Q5"
vg pack -x "$LXGL" -g "$GAM" -Q 5 -t "$T" -o "$PACK"
echo "[$(date)] pack done"

echo "[$(date)] STEP call -v (panel SVs)"
vg call "$LXGL" -k "$PACK" -v "$SVVCF" -t "$T" > "$OUT"
echo "[$(date)] DONE: $(grep -vc '^#' "$OUT") variant lines -> $OUT"
rm -rf "$TMP"
