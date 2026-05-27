#!/bin/bash
#SBATCH --job-name=g0_n231_h
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --array=0-1
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/g0_n231_h_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/g0_n231_h_%A_%a.err

set -euo pipefail
ROOT=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
SIM=$ROOT/control_p80/sims/cov10_n231_g0_s42_hotspots_p80_chr1
OUT_DIR=$ROOT/scratch/g0_n231_h_test
mkdir -p $OUT_DIR

CNS=(cn_full_231_v3qc_v3_filt2
     cn_full_231_v3qc_v3_subsampMedian_refilt2)
CN=${CNS[$SLURM_ARRAY_TASK_ID]}
[[ $CN == *subsamp* ]] && TAG=subsamp || TAG=filt2

echo "=== TASK $SLURM_ARRAY_TASK_ID: cn=$CN  tag=$TAG  sim=cov10_n231_g0 ==="
$PY -u $ROOT/poolfreq/src/sweep_shape_norm_h_only.py \
    --cn-prefix $ROOT/poolfreq/data/$CN/cn_Chr1 \
    --reads $SIM/reads/r1.fq $SIM/reads/r2.fq \
    --sample cov10_n231_g0 \
    --out-prefix $OUT_DIR/${TAG}_cov10_n231_g0 \
    --alphas 0 --threads 8 \
    --counts-cache $OUT_DIR/${TAG}_cov10_n231_g0.counts.npy

echo "DONE $(date)"
