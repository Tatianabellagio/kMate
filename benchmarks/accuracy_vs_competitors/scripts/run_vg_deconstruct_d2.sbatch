#!/bin/bash
#SBATCH --job-name=vg_deconstruct
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/vgdecon_%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/vgdecon_%x_%j.out
set -eo pipefail
# Option 2 (pack-coverage SV AF): deconstruct the SAME d2 graph the pool's pack was
# built on, so the AT (ref+alt node traversals) are in the pack's node-ID space.
# Bridge each record to the panel svidx by REF/ALT SEQUENCE (version-independent).
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vg
G=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.d2.gbz
OUT=benchmarks/accuracy_vs_competitors/work/d2_deconstruct_Chr1.vcf.gz
echo "[$(date)] vg deconstruct d2 (AT traversals, TAIR10#0#Chr1) on $(hostname)"
vg deconstruct -P TAIR10#0#Chr1 -a -e -t "${SLURM_CPUS_PER_TASK:-8}" "$G" | bgzip > "$OUT"
tabix -f -p vcf "$OUT" || true
echo "[$(date)] done: $(zcat "$OUT" | grep -vc '^#') records -> $OUT"
echo "AT populated?"; zcat "$OUT" | grep -v '^#' | head -1 | awk -F'\t' '{print substr($8,1,120)}'