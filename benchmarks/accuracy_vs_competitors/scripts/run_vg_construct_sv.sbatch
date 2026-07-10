#!/bin/bash
#SBATCH --job-name=vg_construct_sv
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/vgconstruct_%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/vgconstruct_%x_%j.out
set -eo pipefail
# -----------------------------------------------------------------------------
# vg-as-SV-AF comparator via the Ash in-silico-pool / vg-issue-4443 method.
#
# Our production graph is minigraph-cactus (.gbz) with NO _alt_ paths, so
# `vg call -v <deconstruct.vcf>` returns 0 variants (even with -z). Fix: build
# a `vg construct -a` graph FROM the panel VCF (which embeds the _alt_ paths
# `vg call -v` needs), map reads to it, then `vg call -v` to get per-SV AD.
#
# Construct from the FULL panel VCF (better read anchoring), but call -v with the
# SV-only VCF (the IDs/positions are a strict subset of the full panel).
# -----------------------------------------------------------------------------
cd /global/scratch/users/tbellg/kmate

POOL=${1:-cov10_n80_g0_s42_hotspots_p80_chr1}
PANEL=${2:-p80}
CHR=${3:-Chr1}
# Which VCF to construct the graph from: "full" (default) or "sv"
CONSTRUCT_FROM=${4:-full}
T=${SLURM_CPUS_PER_TASK:-16}

source ~/miniforge3/etc/profile.d/conda.sh && conda activate vg

REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
FULLVCF=benchmarks/$PANEL/data/pangenome_${PANEL}_chr1.vcf.gz
SVVCF=benchmarks/accuracy_vs_competitors/work/panel_sv_${PANEL}_${CHR}.vcf.gz
READDIR=benchmarks/$PANEL/sims/$POOL/reads

W=benchmarks/accuracy_vs_competitors/work
PFX=$W/c_${POOL}            # index prefix (vg map -d expects PFX.xg + PFX.gcsa)
GAM=$W/${POOL}_construct.gam
PACK=$W/${POOL}_construct.pack
OUT=$W/${POOL}_vgconstruct_sv.vcf

if [ "$CONSTRUCT_FROM" = "sv" ]; then
  CVCF=$SVVCF
else
  CVCF=$FULLVCF
fi

echo "[$(date)] host=$(hostname) cpus=$T construct_from=$CONSTRUCT_FROM"
echo "  REF=$REF"
echo "  CVCF (construct)=$CVCF"
echo "  SVVCF (call -v) =$SVVCF"
echo "  reads=$READDIR/{r1,r2}.fq"

# 1. construct -a (embed _alt_ paths), restrict to CHR ------------------------
echo "[$(date)] STEP1 vg construct -a"
vg construct -a -f -r "$REF" -v "$CVCF" -R "$CHR" -t "$T" > "$PFX.vg"
echo "[$(date)] STEP1 done: $(ls -la $PFX.vg | awk '{print $5}') bytes"

# 2. xg index (for pack + call) ----------------------------------------------
echo "[$(date)] STEP2 vg index -x (xg)"
vg index -x "$PFX.xg" "$PFX.vg" -t "$T"
echo "[$(date)] STEP2 done"

# 3. prune + gcsa index (for mapping) ----------------------------------------
echo "[$(date)] STEP3 vg prune + vg index -g (gcsa)"
vg prune -t "$T" "$PFX.vg" > "$PFX.pruned.vg"
vg index -g "$PFX.gcsa" -k 16 -t "$T" "$PFX.pruned.vg"
echo "[$(date)] STEP3 done"

# 4. map reads (-d uses PFX.xg + PFX.gcsa) -----------------------------------
echo "[$(date)] STEP4 vg map"
vg map -d "$PFX" -f "$READDIR/r1.fq" -f "$READDIR/r2.fq" -t "$T" > "$GAM"
echo "[$(date)] STEP4 done: $(ls -la $GAM | awk '{print $5}') bytes"

# 5. pack (Q5, no vg filter -- vg #4443: filter corrupts AD) ------------------
echo "[$(date)] STEP5 vg pack -Q5"
vg pack -x "$PFX.xg" -g "$GAM" -Q 5 -t "$T" -o "$PACK"
echo "[$(date)] STEP5 done: $(ls -la $PACK | awk '{print $5}') bytes"

# 6. call -v on the SV-only panel VCF ----------------------------------------
echo "[$(date)] STEP6 vg call -v (SV panel)"
vg call "$PFX.xg" -k "$PACK" -v "$SVVCF" -t "$T" > "$OUT"
NV=$(grep -vc '^#' "$OUT" || true)
echo "[$(date)] STEP6 done: $NV variant lines -> $OUT"
echo "=== first 3 records (cols 1-5,9,10) ==="
grep -v '^#' "$OUT" | head -3 | cut -f1-5,9,10 || true
