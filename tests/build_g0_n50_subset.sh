#!/bin/bash
#SBATCH --job-name=g0_n50_build
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/g0_n50_build_%j.out
#SBATCH --error=logs/g0_n50_build_%j.err
mkdir -p logs
set -euo pipefail

# g0 SUBSET pool: pick 50 of the 231 founders (seeded), equal reads each.
# Truth h = 1/50 on the chosen 50, 0 on the other 181. A sparse-h test is more
# discriminating of cn_full identifiability than the perfect 1/231 mix: the EM
# must put mass on the RIGHT 50 and zero elsewhere.

ROOT=/global/scratch/users/tbellg/kmate
WGSIM=/global/home/users/tbellg/miniforge3/envs/pang/bin/wgsim
SAMTOOLS=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/samtools
FOUNDER_DIR=$ROOT/sims/visor_freqk/founder_fastas_231_v3
SIM_DIR=$ROOT/sims/visor_freqk/g0_n50/cov10_g0_n50_chr1
PERFOUNDER=$SIM_DIR/per_founder
mkdir -p $SIM_DIR/reads $PERFOUNDER

# 10x over Chr1 (~30 Mb) / 50 founders = ~20000 pairs/founder (300 bp/pair)
N_PAIRS=20000
SEED_PICK=42

# all founder ids
ALL=$(ls $FOUNDER_DIR/*.chr.fa | xargs -n1 basename | sed 's/.chr.fa//')
# deterministic shuffle -> first 50
CHOSEN=$(echo "$ALL" | sort | awk -v s=$SEED_PICK 'BEGIN{srand(s)} {print rand()"\t"$0}' | sort -k1,1n | cut -f2 | head -50)
echo "$CHOSEN" > $SIM_DIR/chosen_founders.txt
NSEL=$(echo "$CHOSEN" | wc -l)
echo "chosen $NSEL founders -> $SIM_DIR/chosen_founders.txt"

build_one() {
    local fid=$1
    local chr1=$PERFOUNDER/${fid}.chr1.fa
    [[ -s $chr1 ]] || { $SAMTOOLS faidx $FOUNDER_DIR/${fid}.chr.fa Chr1 > $chr1; $SAMTOOLS faidx $chr1; }
    [[ -s $PERFOUNDER/${fid}_1.fq ]] || \
        $WGSIM -N $N_PAIRS -1 150 -2 150 -d 500 -s 50 -e 0.005 -r 0 -R 0 -X 0 \
               -S $(( fid % 100000 + 1 )) \
               $chr1 $PERFOUNDER/${fid}_1.fq $PERFOUNDER/${fid}_2.fq >/dev/null 2>&1
}
export -f build_one; export FOUNDER_DIR PERFOUNDER SAMTOOLS WGSIM N_PAIRS
echo "$CHOSEN" | xargs -n1 -P 16 -I{} bash -c 'build_one "$@"' _ {}

cat $PERFOUNDER/*_1.fq > $SIM_DIR/reads/r1.fq
cat $PERFOUNDER/*_2.fq > $SIM_DIR/reads/r2.fq
echo "pool pairs: $(( $(wc -l < $SIM_DIR/reads/r1.fq) / 4 ))"

# truth h over ALL 231 founders: 1/50 for chosen, 0 else
echo -e "founder\th_truth" > $SIM_DIR/h_truth.tsv
for f in $ALL; do
    if echo "$CHOSEN" | grep -qx "$f"; then v=$(awk -v n=$NSEL 'BEGIN{printf "%.8f",1.0/n}'); else v=0; fi
    echo -e "${f}\t${v}"
done >> $SIM_DIR/h_truth.tsv
echo "wrote $SIM_DIR/h_truth.tsv"

rm -rf $PERFOUNDER
echo "DONE $(date)"
