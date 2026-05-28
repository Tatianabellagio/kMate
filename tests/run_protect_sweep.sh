#!/bin/bash
#SBATCH --job-name=protect_sweep
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --array=0-11
#SBATCH --requeue
#SBATCH --output=logs/protect_sweep_%A_%a.out
#SBATCH --error=logs/protect_sweep_%A_%a.err
mkdir -p logs
set -euo pipefail
ROOT=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
export PYTHONPATH=$ROOT/src:${PYTHONPATH:-}
OUT_DIR=$ROOT/scratch/g0_sweep_h_test
SIM_BASE=$ROOT/sims/visor_freqk/g0_sweep

SIMS=(g0_n10_rep0_cact g0_n50_rep0_cact g0_n50_rep1_bal g0_n50_rep2_pg g0_n200_rep0_rand g0_n231_rep0_rand)
CNS=(cn_full_231_v3qc_v3_subsampProtect1_refilt2 cn_full_231_v3qc_v3_subsampProtect2_refilt2)
TAGS=(protect1 protect2)

N=${#SIMS[@]}
SIM_IDX=$(( SLURM_ARRAY_TASK_ID % N ))
CN_IDX=$(( SLURM_ARRAY_TASK_ID / N ))
SIM=${SIMS[$SIM_IDX]}; CN=${CNS[$CN_IDX]}; TAG=${TAGS[$CN_IDX]}
SIM_DIR=$SIM_BASE/$SIM

echo "[$(date)] EM $TAG on $SIM"
$PY -u $ROOT/src/archive/sweep_shape_norm_h_only.py \
    --cn-prefix $ROOT/data/$CN/cn_Chr1 \
    --reads $SIM_DIR/reads/r1.fq $SIM_DIR/reads/r2.fq \
    --sample $SIM --out-prefix $OUT_DIR/${TAG}_${SIM} \
    --alphas 0 --threads 8 --counts-cache $OUT_DIR/${TAG}_${SIM}.counts.npy
echo "DONE $(date)"
