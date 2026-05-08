#!/bin/bash
#SBATCH --job-name=rgz_bench
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/home/tbellagio/scratch/hapfire_sv/poolfreq/bench_kmc/rgz_bench_%j.out
#SBATCH --error=/home/tbellagio/scratch/hapfire_sv/poolfreq/bench_kmc/rgz_bench_%j.err

# Compare zcat vs rapidgzip as the read-decompression frontend to jellyfish.
# Same inputs as the KMC bench. Goal: see if rapidgzip removes the gzip
# bottleneck (jellyfish was at 555% CPU in the prior bench; if rapidgzip
# unblocks decompression we should see closer to 800% and proportional
# wall-time speedup).
set -uo pipefail
set -x

BENCH=/home/tbellagio/scratch/hapfire_sv/poolfreq/bench_kmc
WORK=${BENCH}/rgz_$SLURM_JOB_ID
mkdir -p $WORK

R1=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz
R2=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz
QUERY_FA=${BENCH}/query_chr1.fa

JELLYFISH=/home/tbellagio/miniforge3/envs/pangenie/bin/jellyfish
RAPIDGZIP=/home/tbellagio/miniforge3/envs/hapfm/bin/rapidgzip
THREADS=8
K=31

# ============================================================================
# CONTROL: zcat | jellyfish (current pipeline) — re-run for clean head-to-head
# ============================================================================
echo "=== ZCAT + JELLYFISH (control) ==="
JF_DB=${WORK}/zcat.jf
JF_OUT=${WORK}/zcat_query.tsv
JF_LOG=${WORK}/zcat.time

/usr/bin/time -v -o $JF_LOG bash -c "
zcat $R1 $R2 | $JELLYFISH count -m $K -s 3G -t $THREADS -C -o $JF_DB /dev/fd/0
$JELLYFISH query $JF_DB -s $QUERY_FA > $JF_OUT
"
echo "ZCAT pipeline:"
grep -E "Elapsed|Maximum resident|Percent of CPU" $JF_LOG

# ============================================================================
# TREATMENT: rapidgzip | jellyfish — same task, just multi-threaded gunzip
# ============================================================================
echo
echo "=== RAPIDGZIP + JELLYFISH ==="
RGZ_DB=${WORK}/rgz.jf
RGZ_OUT=${WORK}/rgz_query.tsv
RGZ_LOG=${WORK}/rgz.time

# rapidgzip can use up to N threads itself; let it have 4 and give jellyfish 4
/usr/bin/time -v -o $RGZ_LOG bash -c "
$RAPIDGZIP -d -P 4 -c $R1 $R2 | $JELLYFISH count -m $K -s 3G -t 4 -C -o $RGZ_DB /dev/fd/0
$JELLYFISH query $RGZ_DB -s $QUERY_FA > $RGZ_OUT
"
echo "RAPIDGZIP pipeline (rapidgzip 4 threads, jellyfish 4 threads):"
grep -E "Elapsed|Maximum resident|Percent of CPU" $RGZ_LOG

# ============================================================================
# TREATMENT 2: rapidgzip 8 threads, single-threaded jellyfish
# (test whether jellyfish or decompression is the binding constraint)
# ============================================================================
echo
echo "=== RAPIDGZIP 8t + JELLYFISH 1t ==="
RGZ8_DB=${WORK}/rgz8.jf
RGZ8_OUT=${WORK}/rgz8_query.tsv
RGZ8_LOG=${WORK}/rgz8.time
/usr/bin/time -v -o $RGZ8_LOG bash -c "
$RAPIDGZIP -d -P 8 -c $R1 $R2 | $JELLYFISH count -m $K -s 3G -t 1 -C -o $RGZ8_DB /dev/fd/0
$JELLYFISH query $RGZ8_DB -s $QUERY_FA > $RGZ8_OUT
"
grep -E "Elapsed|Maximum resident|Percent of CPU" $RGZ8_LOG

# ============================================================================
# TREATMENT 3: rapidgzip 8 threads + jellyfish 8 threads (over-subscribe)
# ============================================================================
echo
echo "=== RAPIDGZIP 8t + JELLYFISH 8t (over-subscribed) ==="
RGZX_DB=${WORK}/rgzx.jf
RGZX_OUT=${WORK}/rgzx_query.tsv
RGZX_LOG=${WORK}/rgzx.time
/usr/bin/time -v -o $RGZX_LOG bash -c "
$RAPIDGZIP -d -P 8 -c $R1 $R2 | $JELLYFISH count -m $K -s 3G -t 8 -C -o $RGZX_DB /dev/fd/0
$JELLYFISH query $RGZX_DB -s $QUERY_FA > $RGZX_OUT
"
grep -E "Elapsed|Maximum resident|Percent of CPU" $RGZX_LOG

# ============================================================================
# AGREEMENT — outputs must match (these are different code paths, must agree)
# ============================================================================
echo
echo "=== AGREEMENT ==="
echo "diff zcat_query vs rgz_query (size, byte-identity):"
ls -la $JF_OUT $RGZ_OUT $RGZ8_OUT $RGZX_OUT
md5sum $JF_OUT $RGZ_OUT $RGZ8_OUT $RGZX_OUT

echo
echo "=== WORKDIR ==="
ls -lh $WORK
