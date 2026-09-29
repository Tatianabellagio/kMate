#!/bin/bash
#SBATCH --job-name=bwa_idx
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=logs/bwa_idx_%j.out
#SBATCH --error=logs/bwa_idx_%j.out

# One-time BWA index for the p231 hapFIRE-vs-kMate speed/accuracy sweep
# (benchmarks/p231 poolsize_depth reads, aligned fresh since those are raw
# VISOR fastq, not the ad-hoc BAM the original single-point benchmark used).
set -euo pipefail
mkdir -p logs
REF=/global/scratch/users/tbellg/pang/ref_xing/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
BWA=/global/home/users/tbellg/miniforge3/envs/kmate/bin/bwa

echo "[$(date)] indexing $REF"
$BWA index "$REF"
echo "[$(date)] done"
ls -la "${REF}".*
