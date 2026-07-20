#!/bin/bash
#SBATCH --job-name=p231_gw_sim
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH --output=logs/gw_sim_%j.out
#SBATCH --error=logs/gw_sim_%j.out

# =============================================================================
# ONE-OFF diagnostic: genome-wide (Chr1-5) reads for the hapFIRE Chr1-only vs
# genome-wide investigation. Stage 1 (mosaics, produces pool_weights.tsv truth)
# + Stage 2 (VISOR reads) only -- Stage 3 (per-variant AF truth) is skipped,
# not needed since this test only scores hapFIRE's ecotype/h accuracy.
# Mirrors benchmarks/p231/scripts/06_run_sim_p231_covarg.sh with CHROMS
# widened to all 5 chromosomes; kept as a separate script (not touching the
# production sim scripts) since this sims dir is throwaway/diagnostic.
#
# Usage: sbatch run_sim_p231_genomewide_test.sh N SEED
# =============================================================================
set -euo pipefail
N_INDIV=${1:?Usage: N SEED}
SEED=${2:?Usage: N SEED}
COVERAGE=10
CHROMS="Chr1 Chr2 Chr3 Chr4 Chr5"

ROOT=/global/scratch/users/tbellg/kmate
CTRL=$ROOT/benchmarks/p231
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
PYTHON=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python

WORK=$ROOT/benchmarks/speed_vs_hapfire/sims_genomewide_test/cov${COVERAGE}_n${N_INDIV}_g0_s${SEED}_p231_genomewide
mkdir -p $WORK logs

CACTUS_DIR=$CTRL/fastas_231
FOUNDERS_META=$ROOT/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1.meta.npz
for f in "$CACTUS_DIR" "$FOUNDERS_META"; do
    [ -e "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

echo "[$(date)] STAGE 1: genome-wide mosaics (--chroms $CHROMS) N=$N_INDIV seed=$SEED"
$PYTHON $ROOT/sims/scripts/make_recomb_mosaics.py \
    --n-indiv $N_INDIV --n-generations 0 --seed $SEED \
    --cactus-dir $CACTUS_DIR --founders-meta $FOUNDERS_META \
    --out-dir $WORK --chroms "$CHROMS" --gen0-no-replace

echo
echo "[$(date)] STAGE 2: VISOR SHORtS at cov ${COVERAGE}x, genome-wide"
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
    if s == 100.0:
        break
    d += 100.0 - s
else:
    raise AssertionError(f'fp fixed-point did not converge; last sum={s!r}')
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
