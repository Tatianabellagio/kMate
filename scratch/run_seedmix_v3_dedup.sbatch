#!/bin/bash
#SBATCH --job-name=sm_v3_dedup
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --array=1-8
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/logs/sm_v3_dedup_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/logs/sm_v3_dedup_%A_%a.err

# cactus_em-v3 global on the 8 SEEDMIX replicates using Trimmomatic+Clumpify-deduped reads
# (NOT BAM-extracted — those lose reads that didn't align to TAIR10).
# Output target: 8 h vectors for cross-replicate variance/bias analysis vs hapFIRE.

set -euo pipefail
PYBIN=/home/tbellagio/miniforge3/envs/hapfm/bin/python
ROOT=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
READS_DIR=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix_trimmed/dedup
OUT_DIR=$ROOT/scratch/seedmix_v3_dedup
mkdir -p $OUT_DIR $ROOT/scratch/logs

S=${SLURM_ARRAY_TASK_ID}
R1=$READS_DIR/seeds-${S}_1_dedup.fq.gz
R2=$READS_DIR/seeds-${S}_2_dedup.fq.gz
OUT=$OUT_DIR/SEEDMIX_S${S}.tsv

echo "=== SEEDMIX_S${S} on $SLURMD_NODENAME, $(date) ==="
echo "R1: $R1"
echo "R2: $R2"
echo "Out: $OUT"
[[ -f $R1 ]] || { echo "MISSING $R1"; exit 1; }
[[ -f $R2 ]] || { echo "MISSING $R2"; exit 1; }

/usr/bin/time -v $PYBIN -u $ROOT/poolfreq/src/per_sample_per_chrom.py \
    --cn-kmer-prefix $ROOT/poolfreq/data/cn_full_231_v3/cn \
    --cn-var $ROOT/poolfreq/data/cn_var_231_v3.cn_var.npz \
    --cn-var-meta $ROOT/poolfreq/data/cn_var_231_v3.meta.npz \
    --reads $R1 $R2 \
    --sample SEEDMIX_S${S} \
    --out $OUT \
    --threads 8

echo "=== SEEDMIX_S${S} done $(date) ==="
ls -lh $OUT $OUT_DIR/SEEDMIX_S${S}.h_per_chrom.npz
