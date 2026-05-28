#!/bin/bash
#SBATCH --job-name=visor_test
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=logs/visor_%j.out
#SBATCH --error=logs/visor_%j.err

# Tier 1 sanity: run our patched per_sample_driver on a few visor_freqk reps
# to verify the new method doesn't regress on the existing benchmark.
#
# Visor_freqk: 231-ecotype uniform pool with a single 1kb deletion injected
# at a known position at known frequency (e.g. rep9 cov50 f50 means deletion
# at Chr1:22373000 in 50% of the 231 ecotypes).
#
# Truth: alt_freq at the deletion position == f/100. Pool h ≈ uniform 1/231.

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate/poolfreq
mkdir -p results/visor_test

PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
DRIVER=src/per_sample_driver.py

# Test 4 combos: rep9 at cov10/20/50 with f50, and rep11 at cov50 f30
for COMBO in \
    "rep9 cov50 f50" \
    "rep9 cov20 f50" \
    "rep9 cov10 f50" \
    "rep11 cov50 f30"
do
    set -- $COMBO
    REP=$1; COV=$2; FREQ=$3
    READS=/global/scratch/users/tbellg/visor_freqk/data/reads_var/del/${REP}/${COV}/var_del_1kb_n231_${FREQ}_err001/all.fq
    OUT=results/visor_test/${REP}_${COV}_${FREQ}_block.tsv
    if [ ! -f "$READS" ]; then
        echo "SKIP $COMBO: $READS missing"; continue
    fi
    if [ -f "$OUT" ]; then
        echo "SKIP $COMBO: output exists"; continue
    fi
    echo ""
    echo "============================================================"
    echo "$COMBO  reads=$(du -h $READS | cut -f1)"
    echo "============================================================"
    $PYTHON -u $DRIVER \
        --cn-kmer-prefix data/cn_full \
        --cn-var data/cn_var_82.cn_var.npz \
        --cn-var-meta data/cn_var_82.meta.npz \
        --reads $READS \
        --sample ${REP}_${COV}_${FREQ} \
        --out $OUT \
        --threads 8 \
        --block-mode window \
        --window-bp 200000
done
echo ""
echo "=== DONE ==="
ls -la results/visor_test/
