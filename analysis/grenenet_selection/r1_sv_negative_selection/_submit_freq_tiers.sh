#!/bin/bash
# Frequency-stratified CDS burden, two tiers (see _sv_burden_freq_bins.py header for why they differ):
#   TIER 1  all 231 founders, CDS split rare (MAC 2-12) / common (MAC 13+).  Private EXCLUDED: it is
#           undetectable in short-read founders (98/151 carry zero private CDS events).
#   TIER 2  the 80 cactus founders only, CDS split private / rare / common. Detection is homogeneous
#           there, so the ultra-rare question can be asked honestly.
# Both on the current main model: pooled SVs, no origin filter, reference TE kept, local-adaptation
# control, plus callrate_sv (the rare-bin burden tracks founder genotyping quality).
#SBATCH --job-name=burden_freq
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-13
#SBATCH --output=logs/burden_freq_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export POOL_SET=sv6 BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export BURDEN_STRATIFY=freq STRATIFY_CTX=CDS
export EXTRA_FIXED="obio1c + obio1c:bio1c + callrate_sv"
# tier 1: 2 bins -> SV3 = sv_CDS_rare, sv_CDS_common, sv_TE, sv_promoter -> 7 tasks (0-6)
# tier 2: 3 bins -> SV3 = +sv_CDS_private                                -> 8 tasks (0-7), cactus only
if [ "$SLURM_ARRAY_TASK_ID" -lt 7 ]; then
  export FREQ_BINS="rare,common"; export SUBSET=""; export TASK=$SLURM_ARRAY_TASK_ID
else
  export FREQ_BINS="private,rare,common"; export SUBSET=cactus; export TASK=$(( SLURM_ARRAY_TASK_ID - 7 ))
fi
echo "[$(date)] host=$(hostname) TIER=$([ "$SLURM_ARRAY_TASK_ID" -lt 7 ] && echo 1 || echo 2) BINS=$FREQ_BINS SUBSET=${SUBSET:-all} TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
