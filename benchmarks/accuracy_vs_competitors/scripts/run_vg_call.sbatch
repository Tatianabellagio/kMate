#!/bin/bash
#SBATCH --job-name=vg_call
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/vgcall_%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/vgcall_%x_%j.out
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
POOL=${1:?pool}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vg
G=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.d2.gbz
PACK=benchmarks/accuracy_vs_competitors/work/${POOL}.pack
OUT=benchmarks/accuracy_vs_competitors/work/${POOL}_vgcall_Chr1.vcf
echo "[$(date)] vg call (all snarls, TAIR10#0#Chr1) on $(hostname): $POOL"
vg call "$G" -k "$PACK" -p TAIR10#0#Chr1 -t "${SLURM_CPUS_PER_TASK:-8}" > "$OUT"
echo "[$(date)] done: $(grep -vc '^#' "$OUT") variant lines -> $OUT"
