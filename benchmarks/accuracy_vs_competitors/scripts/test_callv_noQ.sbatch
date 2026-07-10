#!/bin/bash
#SBATCH --job-name=callv_noQ
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=1:30:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/%x_%j.out
set -eo pipefail
# Isolate the -Q5 hypothesis: re-pack the EXISTING construct GAM with -e (NO -Q, the
# visor recipe) and re-call -v on the -L alts xg. Reuses graph+GAM; no remap/rebuild.
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vg
T=${SLURM_CPUS_PER_TASK:-16}
W=benchmarks/accuracy_vs_competitors/work
POOL=cov10_n80_g0_s42_hotspots_p80_chr1
XGL=$W/c_${POOL}.alts.xg
GAM=$W/${POOL}_construct.gam
SVVCF=$W/panel_sv_p80_Chr1.vcf.gz
PACK=$W/${POOL}_construct_noQ.pack
OUT=$W/${POOL}_vgconstruct_sv_noQ.vcf
echo "[$(date)] host=$(hostname) re-pack -e (NO -Q)"
vg pack -x "$XGL" -g "$GAM" -o "$PACK" -e -t "$T"
echo "[$(date)] vg call -v on no-Q pack"
vg call "$XGL" -k "$PACK" -v "$SVVCF" -t "$T" > "$OUT"
echo "[$(date)] DONE: $(grep -vc '^#' "$OUT") variant lines -> $OUT"
echo "=== alt-depth>0 among SV records ==="
grep -v '^#' "$OUT" | awk -F'\t' '{n=split($9,k,":");for(i=1;i<=n;i++)if(k[i]=="AD"){split($10,v,":");split(v[i],ad,",");a=0;for(j=2;j<=length(ad);j++)a+=ad[j];if(a>0)c++}}END{print c}'
