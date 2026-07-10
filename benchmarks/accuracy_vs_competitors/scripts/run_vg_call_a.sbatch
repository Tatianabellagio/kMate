#!/bin/bash
#SBATCH --job-name=vg_call_a
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/vgcalla_%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/vgcalla_%x_%j.out
set -eo pipefail
# freqk-recipe vg AF: vg call -a (genotype ALL snarls, NOT -v). Emits AD per snarl
# = per panel variant. Reuses the existing pack. Match to panel by (chrom,pos,ref,alt).
cd /global/scratch/users/tbellg/kmate
POOL=${1:?pool}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vg
G=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.d2.gbz
PACK=benchmarks/accuracy_vs_competitors/work/${POOL}.pack
OUT=benchmarks/accuracy_vs_competitors/work/${POOL}_vgcalla_Chr1.vcf
test -f "$PACK" || { echo "missing pack: $PACK"; exit 1; }
echo "[$(date)] vg call -a (ALL snarls) on $(hostname): $POOL"
vg call "$G" -k "$PACK" -a -p TAIR10#0#Chr1 -t "${SLURM_CPUS_PER_TASK:-8}" > "$OUT"
echo "[$(date)] done: $(grep -vc '^#' "$OUT") variant lines -> $OUT"
