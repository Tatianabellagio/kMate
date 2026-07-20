#!/bin/bash
#SBATCH --job-name=gren_sim
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/gren_sim_%j.out
#SBATCH --error=logs/gren_sim_%j.out

# Simulate a pool from the greneNet-derived founder FASTAs (SNP-only genomes),
# for the FAIR hapFIRE benchmark. Reuses the arch3 founders-meta + same seed so
# the founder draw (pool_weights) is IDENTICAL to the arch3 (kMate) p231 sim of
# the same (N,cov,seed) -> matched pools. Only Stage 1 (mosaics -> pool_weights)
# + Stage 2 (VISOR reads) -- no per-variant truth here (scored later against the
# greneNet VCF + pool_weights).
#
# Usage: sbatch run_sim_greneNet.sh N COV SEED
set -euo pipefail
mkdir -p logs
N_INDIV=${1:?Usage: N COV SEED}
COVERAGE=${2:?Usage: N COV SEED}
SEED=${3:?Usage: N COV SEED}
CHROMS="Chr1"

ROOT=/global/scratch/users/tbellg/kmate
BASE=$ROOT/benchmarks/speed_vs_hapfire
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
PYTHON=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python

# arch3 founders-meta (SAME as the kMate p231 sims) -> identical founder draw
FOUNDERS_META=$ROOT/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1.meta.npz
CACTUS_DIR=$BASE/greneNet_fastas
for f in "$CACTUS_DIR" "$FOUNDERS_META"; do
    [ -e "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done
[ "$(ls $CACTUS_DIR/*.chr.fa 2>/dev/null | wc -l)" -eq 231 ] || { echo "ERROR: expected 231 greneNet fastas" >&2; exit 1; }

WORK=$BASE/sims_greneNet/cov${COVERAGE}_n${N_INDIV}_g0_s${SEED}_greneNet_chr1
mkdir -p "$WORK"

echo "[$(date)] STAGE 1: mosaics from greneNet fastas (N=$N_INDIV cov=$COVERAGE seed=$SEED)"
$PYTHON $ROOT/sims/scripts/make_recomb_mosaics.py \
    --n-indiv $N_INDIV --n-generations 0 --seed $SEED \
    --cactus-dir $CACTUS_DIR --founders-meta $FOUNDERS_META \
    --out-dir $WORK --chroms "$CHROMS" --gen0-no-replace

echo
echo "[$(date)] STAGE 2: VISOR SHORtS at cov ${COVERAGE}x"
READS_DIR=$WORK/reads
if [ ! -s ${READS_DIR}/r1.fq ] || [ ! -s ${READS_DIR}/r2.fq ]; then
    REGION_BED=$WORK/region.bed
    declare -A CHROM_LEN=( [Chr1]=30427671 [Chr2]=19698289 [Chr3]=23459830 [Chr4]=18585056 [Chr5]=26975502 )
    > $REGION_BED
    for c in $CHROMS; do
        printf "%s\t1\t%d\t100.0\t100.0\n" "$c" "${CHROM_LEN[$c]}" >> $REGION_BED
    done
    read OTHER_FRAC FIRST_FRAC <<<"$(python3 -c "
N = $N_INDIV
o = round(100.0 / N * 100) / 100.0
d = 100.0 - (N - 1) * o
for _trial in range(50):
    s = sum([d] + [o] * (N - 1))
    if s == 100.0: break
    d += 100.0 - s
else: raise AssertionError('fp fixed-point did not converge')
print(repr(o), repr(d))
")"
    CLONE_DIRS=(); FRACS=(); HAP_DIRS=( $WORK/haps/s_ind* )
    for k in "${!HAP_DIRS[@]}"; do
        CLONE_DIRS+=("${HAP_DIRS[$k]}")
        if [ "$k" = "0" ]; then FRACS+=("$FIRST_FRAC"); else FRACS+=("$OTHER_FRAC"); fi
    done
    echo "  ${#CLONE_DIRS[@]} mosaic clones, others=${OTHER_FRAC}% each, ind001=${FIRST_FRAC}%"
    source /global/home/users/tbellg/miniforge3/etc/profile.d/conda.sh && conda activate kmate
    rm -rf $READS_DIR; mkdir -p $READS_DIR
    VISOR SHORtS \
        -g $REF -s "${CLONE_DIRS[@]}" -b $REGION_BED -o $READS_DIR \
        --coverage $COVERAGE --clonefraction "${FRACS[@]}" \
        --error 0.001 --length 150 --fastq --threads 8 || true
    if [ ! -s ${READS_DIR}/r1.fq ] || [ ! -s ${READS_DIR}/r2.fq ]; then
        echo "ERROR: VISOR SHORtS failed" >&2; exit 1
    fi
fi
echo "[$(date)] DONE -- $WORK"
ls -lh $WORK/reads/r1.fq $WORK/reads/r2.fq $WORK/pool_weights.tsv