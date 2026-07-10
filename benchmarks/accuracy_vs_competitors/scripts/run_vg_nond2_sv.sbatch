#!/bin/bash
#SBATCH --job-name=vg_nond2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/vgnond2_%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/vgnond2_%x_%j.out
set -eo pipefail
# vg-SV done RIGHT: run vg on the NON-d2 graph (the one the panel was deconstructed
# from) so snarl node IDs match the panel/catalog. Then match by traversal node-path.
cd /global/scratch/users/tbellg/kmate
POOL=${1:?pool}; PANEL=${2:-p80}; T=${SLURM_CPUS_PER_TASK:-16}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vg
GDIR=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output
G=$GDIR/pang_1001gplus_82acc.gbz            # NON-d2 (panel node space)
DIST=$GDIR/pang_1001gplus_82acc.dist
MIN=$GDIR/pang_1001gplus_82acc.shortread.min
SNARLS=$GDIR/pang_1001gplus_82acc.snarls
W=benchmarks/accuracy_vs_competitors/work
R1=benchmarks/$PANEL/sims/$POOL/reads/r1.fq; R2=benchmarks/$PANEL/sims/$POOL/reads/r2.fq

# 1. build giraffe indexes for the non-d2 graph (once)
[ -s "$DIST" ] || { echo "[$(date)] build dist"; vg index -j "$DIST" -t "$T" "$G"; }
[ -s "$MIN" ]  || { echo "[$(date)] build minimizer"; vg minimizer -d "$DIST" -o "$MIN" -t "$T" "$G"; }

GAF=$W/${POOL}_nond2.gaf; PACK=$W/${POOL}_nond2.pack; VCF=$W/${POOL}_nond2_vgcalla.vcf
echo "[$(date)] giraffe (non-d2, --rescue-attempts 0: unclipped graph triggers the rescue pathology; AF-equivalent + ~35x faster) on $(hostname): $POOL"
vg giraffe -Z "$G" -d "$DIST" -m "$MIN" --rescue-attempts 0 -f "$R1" -f "$R2" -t "$T" --output-format GAF > "$GAF"
echo "[$(date)] pack"; vg pack -x "$G" -a "$GAF" -Q 5 -o "$PACK" -t "$T"
echo "[$(date)] call -a"; vg call "$G" -k "$PACK" -r "$SNARLS" -a -t "$T" > "$VCF"
rm -f "$GAF" "$PACK"
echo "[$(date)] DONE: $(grep -vc '^#' "$VCF") lines -> $VCF"
