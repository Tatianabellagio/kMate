#!/bin/bash
#SBATCH --job-name=sm_dedup
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --array=1-8
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/arch3/chr1/B1_dedup_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/arch3/chr1/B1_dedup_%A_%a.err
set -euo pipefail

# Re-process SEEDMIX_S{1..8}: add clumpify dedup to xwu's already-pool-seq-trimmed FASTQs.
#
# Why this path (not raw → trim → clumpify): the raw FASTQs at
# /Carnegie/DPB/.../ena_seeds/... are NOT accessible from compute nodes (Carnegie-DPB
# share). xwu's re-trimmed FASTQs use the exact pool-seq Trimmomatic config we want
# (ILLUMINACLIP + SLIDINGWINDOW:4:20 + LEADING:5 + TRAILING:5 + MINLEN:36, from
# xwu/GrENE_net/seed_mix/re-trimmed/commands.sh) and are accessible. So: skip the
# trim step (already done by xwu), just add clumpify dedup (matches preprocess_one.sh
# for founder consistency).
#
# Input  (trimmed-only, no dedup): /home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S{N}-1.{1,2}_P.fq.gz
#   (verified byte-identical to xwu's re-trimmed/S{N}-1.{1,2}_P.fq.gz via md5)
# Output (trim + clumpify dedup):  /home/tbellagio/scratch/pang/grenenet_reads/seed_mix_trimdedup/SEEDMIX_S{N}_{1,2}.dedup.fq.gz
#
# Expected ~35% read reduction post-clumpify (measured on seeds-1: 2.9 GB → 1.9 GB
# in the existing seed_mix_trimmed/dedup/ samples).

S=${SLURM_ARRAY_TASK_ID:?array task id required}

IN_DIR=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix
IN_R1=$IN_DIR/S${S}-1.1_P.fq.gz
IN_R2=$IN_DIR/S${S}-1.2_P.fq.gz

OUT_DIR=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix_trimdedup
mkdir -p $OUT_DIR
OUT_R1=$OUT_DIR/SEEDMIX_S${S}_1.dedup.fq.gz
OUT_R2=$OUT_DIR/SEEDMIX_S${S}_2.dedup.fq.gz

# Idempotent: skip if already done
if [ -s "$OUT_R1" ] && [ -s "$OUT_R2" ]; then
    echo "[$(date)] SEEDMIX_S${S} already dedup'd, skipping"
    ls -lh $OUT_R1 $OUT_R2
    exit 0
fi

[ -s "$IN_R1" ] || { echo "ERROR: missing $IN_R1"; exit 1; }
[ -s "$IN_R2" ] || { echo "ERROR: missing $IN_R2"; exit 1; }

CLUMPIFY=/home/tbellagio/miniforge3/envs/pang/bin/clumpify.sh

# Pre-dedup read counts
echo "[$(date)] === SEEDMIX_S${S}: pre-dedup ==="
PRE_R1=$(zcat $IN_R1 | wc -l | awk '{print $1/4}')
echo "  input R1 reads: $PRE_R1"
echo "  input R1 size:  $(du -h $IN_R1 | cut -f1)"

echo
echo "[$(date)] === SEEDMIX_S${S}: clumpify dedup ==="
# Same flags as preprocess_one.sh: dedupe=t (mark dups), dupesubs=0 (exact match),
# optical=f (no optical-dup detection; we want all PCR dups, not just optical).
$CLUMPIFY in=$IN_R1 in2=$IN_R2 \
    out=$OUT_R1 out2=$OUT_R2 \
    dedupe=t dupesubs=0 optical=f

echo
echo "[$(date)] === SEEDMIX_S${S}: post-dedup ==="
POST_R1=$(zcat $OUT_R1 | wc -l | awk '{print $1/4}')
echo "  output R1 reads: $POST_R1"
echo "  output R1 size:  $(du -h $OUT_R1 | cut -f1)"
echo "  dedup rate:      $(awk "BEGIN{printf \"%.2f\", 100*(1 - $POST_R1/$PRE_R1)}")%"
ls -lh $OUT_R1 $OUT_R2
echo "[$(date)] DONE SEEDMIX_S${S}"
