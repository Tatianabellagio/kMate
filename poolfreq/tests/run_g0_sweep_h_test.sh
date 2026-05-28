#!/bin/bash
#SBATCH --job-name=g0_sweep
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --array=0-17
#SBATCH --requeue
#SBATCH --output=logs/g0_sweep_%A_%a.out
#SBATCH --error=logs/g0_sweep_%A_%a.err

mkdir -p logs
set -euo pipefail
ROOT=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
SIM_BASE=$ROOT/sims/visor_freqk/g0_sweep
OUT_DIR=$ROOT/scratch/g0_sweep_h_test
mkdir -p $SIM_BASE $OUT_DIR

# 9 sim configs:
#   n=10:  rep0=cact-heavy (8C+2PG), rep1=balanced (~3C+7PG), rep2=pg-heavy (1C+9PG)
#   n=50:  rep0=cact-heavy (40C+10PG), rep1=balanced (~17C+33PG), rep2=pg-heavy (5C+45PG)
#   n=200: rep0=random, rep1=random
#   n=231: rep0=all
# Each x 2 cn = 18 SLURM tasks (0..17)

# Pack: [n, rep, n_cactus_or_-1, cn_idx]
CFG=(
  "10  0   8 0" "10  0   8 1"
  "10  1   3 0" "10  1   3 1"
  "10  2   1 0" "10  2   1 1"
  "50  0  40 0" "50  0  40 1"
  "50  1  17 0" "50  1  17 1"
  "50  2   5 0" "50  2   5 1"
  "200 0  -1 0" "200 0  -1 1"
  "200 1  -1 0" "200 1  -1 1"
  "231 0  -1 0" "231 0  -1 1"
)
IDX=$SLURM_ARRAY_TASK_ID
read N REP NC CN_IDX <<< "${CFG[$IDX]}"

CNS=(cn_full_231_v3qc_v3_filt2
     cn_full_231_v3qc_v3_subsampMedian_refilt2)
CN=${CNS[$CN_IDX]}
[[ $CN == *subsamp* ]] && TAG=subsamp || TAG=filt2

# Sim name (matches builder output naming)
if [[ $NC -lt 0 ]]; then
  COMP=rand
elif (( NC*10 >= N*7 )); then
  COMP=cact
elif (( NC*10 <= N*2 )); then
  COMP=pg
else
  COMP=bal
fi
SIM_NAME=g0_n${N}_rep${REP}_${COMP}
SIM_DIR=$SIM_BASE/$SIM_NAME

# Build sim if not already built (only first task per sim does this; rely on r1.fq existence)
if [[ ! -s $SIM_DIR/reads/r1.fq ]]; then
    echo "[$(date)] === Build sim $SIM_NAME (n=$N rep=$REP n_cactus=$NC) ==="
    NC_ARG=""
    [[ $NC -ge 0 ]] && NC_ARG="--n-cactus $NC"
    # Lock with a tmp file to avoid race; if another task is building, wait for r1.fq
    LOCK=$SIM_DIR.building.lock
    mkdir -p $SIM_DIR
    if mkdir $LOCK 2>/dev/null; then
        $PY -u $ROOT/poolfreq/src/build_g0_uniform_sim.py \
            --n $N --rep $REP $NC_ARG \
            --out-base $SIM_BASE
        rmdir $LOCK
    else
        echo "[$(date)] sim build in progress by another task; waiting..."
        while [[ ! -s $SIM_DIR/reads/r1.fq ]] && [[ -d $LOCK ]]; do sleep 30; done
    fi
fi
[[ -s $SIM_DIR/reads/r1.fq ]] || { echo "ERR: $SIM_DIR/reads/r1.fq not built"; exit 1; }

echo "[$(date)] === EM $TAG on $SIM_NAME ==="
$PY -u $ROOT/poolfreq/src/sweep_shape_norm_h_only.py \
    --cn-prefix $ROOT/poolfreq/data/$CN/cn_Chr1 \
    --reads $SIM_DIR/reads/r1.fq $SIM_DIR/reads/r2.fq \
    --sample $SIM_NAME \
    --out-prefix $OUT_DIR/${TAG}_${SIM_NAME} \
    --alphas 0 --threads 8 \
    --counts-cache $OUT_DIR/${TAG}_${SIM_NAME}.counts.npy

echo "DONE $(date)"
