#!/bin/bash
#SBATCH --job-name=compet_km_ld
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/benchmarks/accuracy_vs_competitors/logs/km_ld_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/benchmarks/accuracy_vs_competitors/logs/km_ld_%j.out
set -euo pipefail
mkdir -p /global/scratch/users/tbellg/kmate/benchmarks/accuracy_vs_competitors/{logs,results}
ROOT=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
POOL=cov10_n80_g0_s42_hotspots_p80_chr1
READS=$ROOT/benchmarks/p80/sims/$POOL/reads
KMER=$ROOT/benchmarks/p80/data/kmer_pa_p80_ours_filt2inv/kmer_pa
VAR=$ROOT/benchmarks/p80/data/var_pa_p80
OUT=$ROOT/benchmarks/accuracy_vs_competitors/results/${POOL}_kmate_ldr01.tsv
echo "[$(date)] competitor kMate arm --unit ld on $POOL"
$PY -u $ROOT/src/per_sample_per_chrom.py \
  --kmer-pa-prefix $KMER \
  --var-pa ${VAR}.var_pa.npz --var-called ${VAR}.var_called.npz --var-meta ${VAR}.meta.npz \
  --reads $READS/r1.fq $READS/r2.fq --sample ${POOL}_kmate_ldr01 --out $OUT \
  --threads 8 --chroms Chr1 --unit ld --ld-r2 0.1 --kmer-weight uniform
echo "[$(date)] DONE — $OUT"
