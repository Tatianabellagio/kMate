#!/bin/bash
#SBATCH --job-name=vg_autoidx
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
set -eo pipefail
# FAITHFUL replication of the visor_freqk vg-giraffe SV-AF recipe, applied to the
# kMate benchmark panel SVs. Differences from the broken run_compet_sv path:
#   graph = vg autoindex --workflow giraffe from the SV-ONLY VCF  (not construct -a from FULL panel)
#   pack  = vg pack -e   (NO -Q quality filter -- that filter caused the ~1.7x underestimate)
#   call  = plain vg call (NO -v)
# Usage: sbatch run_vg_sv_autoidx.sbatch <POOL> [PANEL=p80]
cd /global/scratch/users/tbellg/kmate
POOL=${1:?pool}; PANEL=${2:-p80}; T=${SLURM_CPUS_PER_TASK:-32}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vg
W=benchmarks/accuracy_vs_competitors/work
REF_FULL=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
SVVCF=$W/panel_sv_${PANEL}_Chr1.vcf.gz
IDXDIR=$W/autoidx_${PANEL}_sv_Chr1
R1=benchmarks/$PANEL/sims/$POOL/reads/r1.fq; R2=benchmarks/$PANEL/sims/$POOL/reads/r2.fq
CALLS=$W/${POOL}_autoidx_calls.vcf
echo "[$(date)] POOL=$POOL host=$(hostname) T=$T"

# 0. Chr1 reference (cached)
REF=$W/TAIR10_Chr1.fa
[ -s "$REF.fai" ] || { echo "[$(date)] extract Chr1 ref"; samtools faidx "$REF_FULL" Chr1 > "$REF"; samtools faidx "$REF"; }

# 1. autoindex giraffe graph FROM THE SV-ONLY VCF (cached across pools)
if [ ! -s "$IDXDIR/graph.giraffe.gbz" ]; then
  echo "[$(date)] vg autoindex --workflow giraffe from SV VCF ($(zcat $SVVCF|grep -vc '^#') records)"
  mkdir -p "$IDXDIR"
  vg autoindex --workflow giraffe -r "$REF" -v "$SVVCF" -p "$IDXDIR/graph" -t "$T"
else echo "[$(date)] reuse autoindex graph $IDXDIR"; fi
GBZ=$IDXDIR/graph.giraffe.gbz; DIST=$IDXDIR/graph.dist
MIN=$(ls "$IDXDIR"/graph*.min 2>/dev/null | head -1)
ZIP=$(ls "$IDXDIR"/graph*.zipcodes 2>/dev/null | head -1 || true)
echo "[$(date)] MIN=$MIN ZIP=${ZIP:-none}"

# 2. stage indexes node-local (the I/O fix), then giraffe
LOC=${TMPDIR:-/tmp}/autoidx; mkdir -p "$LOC"
cp "$GBZ" "$DIST" "$MIN" $ZIP "$LOC/" 2>/dev/null
LGBZ=$LOC/$(basename "$GBZ"); LDIST=$LOC/$(basename "$DIST"); LMIN=$LOC/$(basename "$MIN")
LZIP=""; [ -n "$ZIP" ] && LZIP="-z $LOC/$(basename "$ZIP")"
GAM=$LOC/${POOL}.gam; PACK=$LOC/${POOL}.pack
echo "[$(date)] vg giraffe (node-local, visor recipe)"
vg giraffe -Z "$LGBZ" -d "$LDIST" -m "$LMIN" $LZIP -f "$R1" -f "$R2" -t "$T" -p > "$GAM"

# 3. pack -e (NO -Q)
echo "[$(date)] vg pack -e (no -Q)"
vg pack -x "$LGBZ" -g "$GAM" -o "$PACK" -e -t "$T"

# 4. plain vg call (NO -v)
echo "[$(date)] vg call (no -v)"
vg call "$LGBZ" -k "$PACK" -t "$T" > "$CALLS"
echo "[$(date)] DONE: $(grep -vc '^#' "$CALLS") call records -> $CALLS"
