#!/bin/bash
#SBATCH --job-name=compet_sv
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/csv_%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/csv_%x_%j.out
set -eo pipefail
# SV-side competitor for ONE p80 pool: vg on the NON-d2 graph (panel node space) ->
# pack -> call -a -> traversal-matched per-panel-SV AF (svidx,est). Idempotent.
# Usage: sbatch run_compet_sv.sbatch <POOL> [PANEL=p80]
#
# SPEED FIX (2026-06-24): vg giraffe mmaps its indexes and random-accesses them per
# read; reading them off Lustre turns every page fault into a network RPC -> ~85x
# slowdown (jobs died at the 24h wall). Fix = stage indexes to NODE-LOCAL $TMPDIR
# first, and use the ZIPCODE minimizer (.withzip.min + .zipcodes) not the plain one.
cd /global/scratch/users/tbellg/kmate
POOL=${1:?pool}; PANEL=${2:-p80}; T=${SLURM_CPUS_PER_TASK:-16}
source ~/miniforge3/etc/profile.d/conda.sh
SC=benchmarks/accuracy_vs_competitors/scripts
W=benchmarks/accuracy_vs_competitors/work
GDIR=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output
PFX=pang_1001gplus_82acc                     # NON-d2 (panel node space)
PANELSV=$W/panel_sv_${PANEL}_Chr1.vcf.gz
R1=benchmarks/$PANEL/sims/$POOL/reads/r1.fq; R2=benchmarks/$PANEL/sims/$POOL/reads/r2.fq
OUT=$W/vg_sv_${POOL}.tsv
VCF=$W/${POOL}_nond2_vgcalla.vcf
echo "[$(date)] POOL=$POOL PANEL=$PANEL host=$(hostname) cpus=$T"

if [ -s "$OUT" ]; then echo "[$(date)] skip (have $OUT)"; exit 0; fi

if [ ! -s "$VCF" ]; then
  conda activate vg
  # --- stage mmap'd indexes to node-local disk (the speed fix) ---
  LOC=${TMPDIR:-/tmp}/giridx_nond2; mkdir -p "$LOC"
  echo "[$(date)] staging non-d2 indexes -> $LOC"; t0=$SECONDS
  cp $GDIR/$PFX.gbz $GDIR/$PFX.dist $GDIR/$PFX.shortread.withzip.min \
     $GDIR/$PFX.shortread.zipcodes $GDIR/$PFX.snarls "$LOC/"
  echo "[$(date)] staged in $((SECONDS-t0))s"
  G=$LOC/$PFX.gbz; DIST=$LOC/$PFX.dist; MIN=$LOC/$PFX.shortread.withzip.min
  ZIP=$LOC/$PFX.shortread.zipcodes; SNARLS=$LOC/$PFX.snarls
  GAF=$LOC/${POOL}_nond2.gaf; PACK=$LOC/${POOL}_nond2.pack
  echo "[$(date)] giraffe (non-d2, zipcode, node-local, --rescue-attempts 0): $POOL"
  vg giraffe -Z "$G" -d "$DIST" -m "$MIN" -z "$ZIP" --rescue-attempts 0 \
     -f "$R1" -f "$R2" -t "$T" -p --output-format GAF > "$GAF"
  echo "[$(date)] pack"; vg pack -x "$G" -a "$GAF" -Q 5 -o "$PACK" -t "$T"
  echo "[$(date)] call -a"; vg call "$G" -k "$PACK" -r "$SNARLS" -a -t "$T" > "$VCF"
  rm -f "$GAF" "$PACK"
else echo "[$(date)] reuse vcf $VCF"; fi

echo "[$(date)] traversal-match -> $OUT"
conda activate kmate
python $SC/vg_nond2_sv_af.py --vg-vcf "$VCF" --panel-sv-vcf "$PANELSV" --out "$OUT" --svlen 50
rm -f "$VCF"     # 1.5GB each; drop after conversion
echo "[$(date)] DONE SV-side $POOL"
