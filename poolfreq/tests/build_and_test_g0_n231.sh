#!/bin/bash
#SBATCH --job-name=g0_n231
#SBATCH --partition=bse
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/g0_n231_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/g0_n231_%j.err

set -euo pipefail
ROOT=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
WGSIM=/home/tbellagio/miniforge3/envs/pang/bin/wgsim
SAMTOOLS=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/samtools

SIM_DIR=$ROOT/sims/visor_freqk/g0_n231/cov10_g0_n231_chr1
FOUNDER_DIR=$ROOT/sims/visor_freqk/founder_fastas_231_v3
PERFOUNDER=$SIM_DIR/per_founder
mkdir -p $SIM_DIR/reads $PERFOUNDER

echo "=== Build g0_n231 sim ($(date)) ==="
# Chr1 = 30,003,408 bp; cov=10× → 300M bases / (2×150 bp/pair) ~ 1M pairs; /231 founders → 4329 pairs each
N_PAIRS=4329

# Per-founder Chr1 extract + wgsim, parallel
build_one() {
    local fid=$1
    local chr1=$PERFOUNDER/${fid}.chr1.fa
    if [[ ! -s $chr1 ]]; then
        $SAMTOOLS faidx $FOUNDER_DIR/${fid}.chr.fa Chr1 > $chr1
        $SAMTOOLS faidx $chr1
    fi
    if [[ ! -s $PERFOUNDER/${fid}_1.fq ]]; then
        $WGSIM -N $N_PAIRS -1 150 -2 150 -d 500 -s 50 -e 0.005 -r 0 -R 0 -X 0 \
               -S $(( fid % 100000 + 1 )) \
               $chr1 $PERFOUNDER/${fid}_1.fq $PERFOUNDER/${fid}_2.fq >/dev/null 2>&1
    fi
}
export -f build_one
export FOUNDER_DIR PERFOUNDER SAMTOOLS WGSIM N_PAIRS

FIDS=$(ls $FOUNDER_DIR/*.chr.fa | xargs -n1 basename | sed 's/.chr.fa//')
N=$(echo "$FIDS" | wc -l)
echo "  building $N founder read sets in parallel..."
echo "$FIDS" | xargs -n1 -P 16 -I{} bash -c 'build_one "$@"' _ {}

echo "  concatenating per-founder reads → pool r1.fq, r2.fq"
cat $PERFOUNDER/*_1.fq > $SIM_DIR/reads/r1.fq
cat $PERFOUNDER/*_2.fq > $SIM_DIR/reads/r2.fq
echo "  total pool reads: $(wc -l < $SIM_DIR/reads/r1.fq | awk '{print $1/4}') pairs"

# Truth h: uniform 1/231 for all
echo -e "founder\th_truth" > $SIM_DIR/h_truth.tsv
for f in $(echo "$FIDS"); do
    echo -e "${f}\t$(awk -v n=$N 'BEGIN{printf "%.6f", 1.0/n}')"
done >> $SIM_DIR/h_truth.tsv
echo "  wrote $SIM_DIR/h_truth.tsv"

# Clean per-founder intermediates to save space
echo "  cleaning per-founder intermediates..."
rm -rf $PERFOUNDER

echo
echo "=== Run cactus_em filt2 baseline ==="
$PY -u $ROOT/poolfreq/src/sweep_shape_norm_h_only.py \
    --cn-prefix $ROOT/poolfreq/data/cn_full_231_v3qc_v3_filt2/cn_Chr1 \
    --reads $SIM_DIR/reads/r1.fq $SIM_DIR/reads/r2.fq \
    --sample g0_n231 \
    --out-prefix $SIM_DIR/filt2_g0_n231 \
    --alphas 0 --threads 8 \
    --counts-cache $SIM_DIR/filt2_g0_n231.counts.npy

echo
echo "=== Run cactus_em subsamp_raw cn ==="
$PY -u $ROOT/poolfreq/src/sweep_shape_norm_h_only.py \
    --cn-prefix $ROOT/poolfreq/data/cn_full_231_v3qc_v3_subsampMedian_refilt2/cn_Chr1 \
    --reads $SIM_DIR/reads/r1.fq $SIM_DIR/reads/r2.fq \
    --sample g0_n231 \
    --out-prefix $SIM_DIR/subsamp_g0_n231 \
    --alphas 0 --threads 8 \
    --counts-cache $SIM_DIR/subsamp_g0_n231.counts.npy

echo "DONE $(date)"
