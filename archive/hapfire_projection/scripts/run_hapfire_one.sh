#!/bin/bash
#SBATCH --job-name=hapfire_sv
#SBATCH --output=/home/tbellagio/scratch/hapfire_sv/logs/hapfire_%x_%j.out
#SBATCH --error=/home/tbellagio/scratch/hapfire_sv/logs/hapfire_%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G

set -euo pipefail

# Args:
#   $1 = BAM path (Chr1-renamed if needed)
#   $2 = output prefix (full path including basename)
BAM="$1"
OUT_PREFIX="$2"

VCF=/home/tbellagio/scratch/hapfire_sv/data/vcf/greneNet_Chr1_only.vcf
REF=/home/tbellagio/scratch/visor_freqk/data/reference/Chr1.fa
HAPFIRE_DIR=/carnegie/nobackup/scratch/xwu/haplotype_frequency_estimation/hapFIRE_sourcecode
HAPFIRE=${HAPFIRE_DIR}/hapFIRE.py
HARP_BIN=${HAPFIRE_DIR}/bin

PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python

# harp must be on PATH
export PATH="${HARP_BIN}:${PATH}"

# Working dir: hapFIRE writes intermediate files alongside output prefix
mkdir -p "$(dirname "${OUT_PREFIX}")"
cd "$(dirname "${OUT_PREFIX}")"

# index BAM if needed
if [[ ! -f "${BAM}.bai" ]]; then
  mamba run -n pang samtools index "${BAM}"
fi

# index reference if needed
if [[ ! -f "${REF}.fai" ]]; then
  mamba run -n pang samtools faidx "${REF}"
fi

echo "[$(date)] running hapFIRE on ${BAM} -> $(basename "${OUT_PREFIX}")"
${PYTHON} ${HAPFIRE} \
    -v "${VCF}" \
    -b "${BAM}" \
    -f "${REF}" \
    -o "$(basename "${OUT_PREFIX}")"

echo "[$(date)] hapFIRE done: ${OUT_PREFIX}"
ls -la "$(dirname "${OUT_PREFIX}")"/$(basename "${OUT_PREFIX}")*
