#!/bin/bash
#SBATCH --job-name=h_acc_ld
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/benchmarks/h_accuracy/logs/h_%A_%a.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/benchmarks/h_accuracy/logs/h_%A_%a.out

# Founder-h recovery under the NEW production estimator, --unit chrom (one h per
# chromosome = the founder-decomposition estimand). Production in-house panel +
# per_founder + uniform + haploblock eps=0. Writes <out>.h_per_chrom.npz for
# score_h_vs_truth.py. (h is a one-h estimand → --unit chrom, not ld.)
#
# Usage: sbatch --array=0-1 run_h_unit_chrom.sbatch
mkdir -p /global/scratch/users/tbellg/kmate/benchmarks/h_accuracy/logs
set -euo pipefail
ROOT=/global/scratch/users/tbellg/kmate
PYTHON=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
DRIVER=$ROOT/src/per_sample_per_chrom.py
KMER=$ROOT/data/kmer_pa_231_arch3_filt2inv/kmer_pa
OUT=$ROOT/benchmarks/h_accuracy

REGIMES=(n231_g0 n50_g0)
SUBDIRS=(cov10_n231_g0_s42_hotspots_p231_chr1 cov10_n50_g0_s42_hotspots_p231_chr1)
i=${SLURM_ARRAY_TASK_ID:-0}
REG=${REGIMES[$i]}; SUB=${SUBDIRS[$i]}
READS=$ROOT/benchmarks/p231/sims/$SUB/reads

echo "[$(date)] h-recovery --unit chrom  regime=$REG"
$PYTHON -u $DRIVER \
    --kmer-pa-prefix $KMER \
    --reads $READS/r1.fq $READS/r2.fq \
    --sample ldr01_${REG} --out $OUT/ldr01_${REG}.tsv \
    --threads 8 --chroms Chr1 \
    --unit chrom --h-only --kmer-weight uniform
echo "[$(date)] DONE — $OUT/ldr01_${REG}.h_per_chrom.npz"
