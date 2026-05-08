#!/bin/bash
#SBATCH --job-name=kmc_bench
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/home/tbellagio/scratch/hapfire_sv/poolfreq/bench_kmc/bench_%j.out
#SBATCH --error=/home/tbellagio/scratch/hapfire_sv/poolfreq/bench_kmc/bench_%j.err

# Benchmark jellyfish vs KMC on SEEDMIX_S1 reads, Chr1 query set (20.8M k-mers).
# Same reads, same query, same threads (8). Compare wall time, peak RSS, and
# output identity (every k-mer's count must agree).
set -uo pipefail
set -x

BENCH=/home/tbellagio/scratch/hapfire_sv/poolfreq/bench_kmc
WORK=${BENCH}/work_$SLURM_JOB_ID
mkdir -p $WORK

R1=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz
R2=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz
QUERY_FA=${BENCH}/query_chr1.fa

JELLYFISH=/home/tbellagio/miniforge3/envs/pangenie/bin/jellyfish
KMC=/home/tbellagio/miniforge3/envs/BIOS424/bin/kmc
KMC_TOOLS=/home/tbellagio/miniforge3/envs/BIOS424/bin/kmc_tools

THREADS=8
K=31

# ============================================================================
# 1. JELLYFISH (current production)
# ============================================================================
echo "=== JELLYFISH ==="
JF_DB=${WORK}/reads.jf
JF_OUT=${WORK}/jf_query.tsv
JF_LOG=${WORK}/jf.time

/usr/bin/time -v -o $JF_LOG bash -c "
zcat $R1 $R2 | $JELLYFISH count -m $K -s 3G -t $THREADS -C -o $JF_DB /dev/fd/0
$JELLYFISH query $JF_DB -s $QUERY_FA > $JF_OUT
"

echo "JELLYFISH time/memory:"
grep -E "Elapsed|Maximum resident" $JF_LOG
ls -lh $JF_DB $JF_OUT

# ============================================================================
# 2. KMC 3
# ============================================================================
echo
echo "=== KMC 3 ==="
KMC_DB=${WORK}/reads_kmc
KMC_QUERY_DB=${WORK}/query_kmc
KMC_INTER=${WORK}/intersect_kmc
KMC_OUT=${WORK}/kmc_query.tsv
KMC_LOG=${WORK}/kmc.time
KMC_TMP=${WORK}/kmc_tmp
mkdir -p $KMC_TMP

# Step 1: count k-mers in reads (-ci0 keeps even k=1 counts)
# kmc reads gz natively via @-list (one path per line) or stdin
READS_LIST=${WORK}/reads.list
echo $R1 > $READS_LIST
echo $R2 >> $READS_LIST

/usr/bin/time -v -o ${KMC_LOG}.count $KMC -k$K -ci1 -t$THREADS -m48 \
    @$READS_LIST $KMC_DB $KMC_TMP

# Step 2: build query k-mer database (counts will be 1 each)
/usr/bin/time -v -o ${KMC_LOG}.query_db $KMC -k$K -ci1 -fa -t$THREADS -m8 \
    $QUERY_FA $KMC_QUERY_DB $KMC_TMP

# Step 3: intersect, taking counts from reads_kmc (-ocfirst)
# (default counter combination is min; -ocfirst takes from first input)
/usr/bin/time -v -o ${KMC_LOG}.intersect $KMC_TOOLS -t$THREADS simple \
    $KMC_DB -ci1 $KMC_QUERY_DB -ci1 \
    intersect $KMC_INTER -ocfirst

# Step 4: dump intersection to text
$KMC_TOOLS transform $KMC_INTER dump $KMC_OUT

echo "KMC steps time/memory:"
for stage in count query_db intersect; do
    echo "-- $stage --"
    grep -E "Elapsed|Maximum resident" ${KMC_LOG}.${stage}
done
ls -lh ${KMC_DB}.* $KMC_OUT

# ============================================================================
# 3. AGREEMENT CHECK
# ============================================================================
echo
echo "=== AGREEMENT CHECK ==="
/home/tbellagio/miniforge3/envs/hapfm/bin/python ${BENCH}/compare_outputs.py \
    --jf $JF_OUT --kmc $KMC_OUT --query $QUERY_FA

echo
echo "=== DONE — workdir: $WORK ==="
ls -lh $WORK
