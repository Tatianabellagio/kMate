#!/bin/bash
#SBATCH --job-name=g0_n231_h
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --array=0-1
#SBATCH --requeue
#SBATCH --output=logs/g0_n231_h_%A_%a.out
#SBATCH --error=logs/g0_n231_h_%A_%a.err

mkdir -p logs
set -euo pipefail
ROOT=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
SIM=$ROOT/benchmarks/p80/sims/cov10_n231_g0_s42_hotspots_p80_chr1
OUT_DIR=$ROOT/scratch/g0_n231_h_test
mkdir -p $OUT_DIR

CNS=(kmer_pa_231_v3qc_v3_filt2
     kmer_pa_231_v3qc_v3_subsampMedian_refilt2)
CN=${CNS[$SLURM_ARRAY_TASK_ID]}
[[ $CN == *subsamp* ]] && TAG=subsamp || TAG=filt2

echo "=== TASK $SLURM_ARRAY_TASK_ID: kmer_pa=$CN  tag=$TAG  sim=cov10_n231_g0 ==="
$PY -u $ROOT/src/sweep_shape_norm_h_only.py \
    --kmer_pa-prefix $ROOT/data/$CN/kmer_pa_Chr1 \
    --reads $SIM/reads/r1.fq $SIM/reads/r2.fq \
    --sample cov10_n231_g0 \
    --out-prefix $OUT_DIR/${TAG}_cov10_n231_g0 \
    --alphas 0 --threads 8 \
    --counts-cache $OUT_DIR/${TAG}_cov10_n231_g0.counts.npy

echo "DONE $(date)"
