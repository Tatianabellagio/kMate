#!/bin/bash
#SBATCH --job-name=rgz_bench2
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=/home/tbellagio/scratch/hapfire_sv/poolfreq/bench_kmc/rgz_bench2_%j.out
#SBATCH --error=/home/tbellagio/scratch/hapfire_sv/poolfreq/bench_kmc/rgz_bench2_%j.err

# Re-bench: rapidgzip processes one positional file, so to handle paired
# FASTQs we either cat the gz files (gzip is concat-safe) or run rapidgzip
# twice. Test both invocations.
set -uo pipefail
set -x

BENCH=/home/tbellagio/scratch/hapfire_sv/poolfreq/bench_kmc
WORK=${BENCH}/rgz2_$SLURM_JOB_ID
mkdir -p $WORK

R1=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz
R2=/home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz
QUERY_FA=${BENCH}/query_chr1.fa

JELLYFISH=/home/tbellagio/miniforge3/envs/pangenie/bin/jellyfish
RAPIDGZIP=/home/tbellagio/miniforge3/envs/hapfm/bin/rapidgzip
THREADS=8

# === A: cat R1 R2 | rapidgzip -d (gzip is concatenation-safe) ===
echo "=== A: cat R1 R2 | rapidgzip 8t | jellyfish 8t ==="
A_LOG=${WORK}/A.time
A_OUT=${WORK}/A_query.tsv
A_DB=${WORK}/A.jf
/usr/bin/time -v -o $A_LOG bash -c "
cat $R1 $R2 | $RAPIDGZIP -d -P 8 -c | $JELLYFISH count -m 31 -s 3G -t 8 -C -o $A_DB /dev/fd/0
$JELLYFISH query $A_DB -s $QUERY_FA > $A_OUT
"
grep -E "Elapsed|Maximum resident|Percent of CPU" $A_LOG

# === B: sequential rapidgzip then concat ===
echo
echo "=== B: rapidgzip R1; rapidgzip R2 | jellyfish 8t (sequential) ==="
B_LOG=${WORK}/B.time
B_OUT=${WORK}/B_query.tsv
B_DB=${WORK}/B.jf
/usr/bin/time -v -o $B_LOG bash -c "
{ $RAPIDGZIP -d -P 8 -c $R1; $RAPIDGZIP -d -P 8 -c $R2; } | $JELLYFISH count -m 31 -s 3G -t 8 -C -o $B_DB /dev/fd/0
$JELLYFISH query $B_DB -s $QUERY_FA > $B_OUT
"
grep -E "Elapsed|Maximum resident|Percent of CPU" $B_LOG

# === COMPARE to known-good zcat output ===
echo
echo "=== AGREEMENT vs zcat (control from prior bench) ==="
md5sum ${BENCH}/rgz_59320/zcat_query.tsv $A_OUT $B_OUT
ls -lh ${BENCH}/rgz_59320/zcat_query.tsv $A_OUT $B_OUT

# Quick value-level check on A
/home/tbellagio/miniforge3/envs/hapfm/bin/python -c "
zcat = {}
with open('${BENCH}/rgz_59320/zcat_query.tsv') as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 2:
            zcat[parts[0]] = int(parts[1])
A = {}
with open('$A_OUT') as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 2:
            A[parts[0]] = int(parts[1])
print(f'zcat sum: {sum(zcat.values()):,}')
print(f'A    sum: {sum(A.values()):,}')
n_diff = sum(1 for k in zcat if zcat.get(k,0) != A.get(k,0))
print(f'kmers with different counts: {n_diff:,}')
"
