#!/bin/bash
#SBATCH --job-name=sim_h
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --array=0-5
#SBATCH --requeue
#SBATCH --output=logs/sim_h_%A_%a.out
#SBATCH --error=logs/sim_h_%A_%a.err

mkdir -p logs
set -euo pipefail
ROOT=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
SIM_BASE=$ROOT/sims/visor_freqk/pool_sweep_82_recomb
OUT_DIR=$ROOT/scratch/sim_h_test
mkdir -p $OUT_DIR

# Job grid: 3 sims x 2 kmer_pa matrices = 6 jobs
SIMS=(cov10_n200_g1_s42_hotspots_p231_chr1
      cov10_n50_g1_s42_hotspots_p231_chr1
      cov10_n50_g3_s42_hotspots_p231_chr1)
CNS=(kmer_pa_231_v3qc_v3_filt2
     kmer_pa_231_v3qc_v3_subsampMedian_refilt2)

SIM_IDX=$(( SLURM_ARRAY_TASK_ID / 2 ))
CN_IDX=$(( SLURM_ARRAY_TASK_ID % 2 ))
SIM=${SIMS[$SIM_IDX]}
CN=${CNS[$CN_IDX]}

if [[ $CN == *subsamp* ]]; then
    TAG=subsamp
else
    TAG=filt2
fi

echo "=== TASK $SLURM_ARRAY_TASK_ID: sim=$SIM  kmer_pa=$CN  tag=$TAG ==="
$PY -u $ROOT/src/sweep_shape_norm_h_only.py \
    --kmer_pa-prefix $ROOT/data/$CN/kmer_pa_Chr1 \
    --reads $SIM_BASE/$SIM/reads/r1.fq $SIM_BASE/$SIM/reads/r2.fq \
    --sample $SIM \
    --out-prefix $OUT_DIR/${TAG}_${SIM} \
    --alphas 0 --threads 8 \
    --counts-cache $OUT_DIR/${TAG}_${SIM}.counts.npy

echo "DONE $(date)"
