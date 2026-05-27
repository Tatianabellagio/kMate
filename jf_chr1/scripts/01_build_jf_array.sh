#!/bin/bash
#SBATCH --job-name=jf_build
#SBATCH --partition=bse
#SBATCH --array=1-82%50
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/jf_chr1/logs/build_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/jf_chr1/logs/build_%A_%a.err

set -euo pipefail
ROOT=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/jf_chr1
JF=/home/tbellagio/miniforge3/envs/pangenie/bin/jellyfish
SAMTOOLS=/home/tbellagio/miniforge3/envs/BIOS424/bin/samtools
CHRDIR=/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only

FID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $ROOT/scripts/asm_ids.txt)
OUT=$ROOT/jf/$FID.jf

# Skip if already built
if [ -s "$OUT" ]; then
  echo "[$(date +%T)] $FID: existing .jf found, skipping"
  exit 0
fi

# Extract Chr1, count canonical k-mers
TMP=$(mktemp -d -p /tmp jf_${FID}_XXXX)
trap "rm -rf $TMP" EXIT
$SAMTOOLS faidx $CHRDIR/$FID.chr.fa Chr1 > $TMP/$FID.chr1.fa
$JF count -m 31 -s 50M -t 2 -C -o $OUT $TMP/$FID.chr1.fa
$JF stats $OUT | head -3
echo "[$(date +%T)] $FID done"
