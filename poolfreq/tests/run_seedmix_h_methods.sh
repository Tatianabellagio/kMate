#!/bin/bash
#SBATCH --job-name=sm_h
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=2:00:00
#SBATCH --array=0-23
#SBATCH --requeue
#SBATCH --output=logs/sm_h_%A_%a.out
#SBATCH --error=logs/sm_h_%A_%a.err
mkdir -p logs
set -euo pipefail
ROOT=/global/scratch/users/tbellg/hapfire_sv
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
export PYTHONPATH=$ROOT/poolfreq/src:${PYTHONPATH:-}
READS=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix_trimdedup
OUT=$ROOT/scratch/seedmix_h_test; mkdir -p $OUT

CNS=(cn_full_231_v3qc_v3_filt2 cn_full_231_v3qc_v3_subsampMedian_refilt2 cn_full_231_v3qc_v3_subsampProtect1_refilt2)
TAGS=(filt2 subsamp protect1)
NMETH=3
REP=$(( SLURM_ARRAY_TASK_ID % 8 + 1 ))     # S1..S8
M=$(( SLURM_ARRAY_TASK_ID / 8 ))           # 0..2
CN=${CNS[$M]}; TAG=${TAGS[$M]}

R1=$READS/SEEDMIX_S${REP}_1.dedup.fq.gz
R2=$READS/SEEDMIX_S${REP}_2.dedup.fq.gz
echo "[$(date)] EM $TAG on SEEDMIX_S$REP"
$PY -u $ROOT/poolfreq/src/archive/sweep_shape_norm_h_only.py \
    --cn-prefix $ROOT/poolfreq/data/$CN/cn_Chr1 \
    --reads $R1 $R2 --sample SEEDMIX_S$REP \
    --out-prefix $OUT/${TAG}_S${REP} \
    --alphas 0 --threads 8 --counts-cache $OUT/${TAG}_S${REP}.counts.npy
echo "DONE $(date)"
