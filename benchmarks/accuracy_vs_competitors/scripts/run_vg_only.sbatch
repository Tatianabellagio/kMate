#!/bin/bash
#SBATCH --job-name=vg_only
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/vgonly_%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/vgonly_%x_%j.out
set -eo pipefail
# Relocatable: run from the submission dir (submit_vg_p231_s42.sh cd's to repo root
# before sbatch, so SLURM_SUBMIT_DIR is the repo root) instead of a hardcoded path.
cd "${SLURM_SUBMIT_DIR:-.}"
POOL=${1:?pool}; PANEL=${2:-p80}; PRESET=${3:-}
# Canonical benchmark config: rescue-OFF (matches run_compet_snp.sbatch). Default
# giraffe is rescue-bound -> ~35x slower + times out on the dense SNP-construct
# graph; rescue-off is AF-equivalent. Override by exporting VG_RESCUE before sbatch.
export VG_RESCUE=${VG_RESCUE:-0}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vg
/usr/bin/time -v bash benchmarks/accuracy_vs_competitors/scripts/run_vg.sh "$POOL" "$PANEL" "${SLURM_CPUS_PER_TASK:-8}" "$PRESET"
