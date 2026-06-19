#!/bin/bash
#SBATCH --job-name=p80_b_sim
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/06_sim_%j.out
#SBATCH --error=logs/06_sim_%j.err

# =============================================================================
# Phase B -- recombinant cov10 sim on the p80 panel (Chr1 only).
#
# Usage:
#   sbatch 06_run_sim_p80.sh N_INDIV N_GEN [SEED=42] [SELFING_RATE=0.0]
# e.g.
#   sbatch 06_run_sim_p80.sh 50 1   # n50_g1
#   sbatch 06_run_sim_p80.sh 50 3   # n50_g3
#   sbatch 06_run_sim_p80.sh 80 1   # n80_g1 (stretch; full-panel pool)
#   sbatch 06_run_sim_p80.sh 50 3 42 0.97   # n50_g3 with 97% selfing (-> _self97 tag)
#
# SELFING_RATE>0 models A. thaliana selfing: each offspring is a clonal copy of
# one parent w.p. SELFING_RATE, else a recombinant outcross. Output dir gets a
# _self<pct> tag so it never collides with the forced-outcross (legacy) sims.
#
# Self-contained: clones of make_recomb_mosaics.py + VISOR SHORtS + truth, all
# pointing at benchmarks/p80/data/* artifacts. No v3 paths.
# =============================================================================
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

N_INDIV=${1:?Usage: sbatch 06_run_sim_p80.sh N_INDIV N_GEN [SEED=42] [SELFING_RATE=0.0]}
N_GEN=${2:?Usage: sbatch 06_run_sim_p80.sh N_INDIV N_GEN [SEED=42] [SELFING_RATE=0.0]}
SEED=${3:-42}
SELFING_RATE=${4:-0.0}
COVERAGE=10  # locked: cov10 for this experiment
CHROMS="Chr1"

CTRL=/global/scratch/users/tbellg/kmate/benchmarks/p80
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
PYTHON=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
SCRIPTS=$CTRL/scripts

# Selfing tag/flag: empty for the legacy (rate 0.0) sims so existing dir paths
# and behaviour are byte-identical; _self<pct> otherwise.
SELF_TAG=""; SELF_FLAG=""
if [ "$(python3 -c "print(float('$SELFING_RATE')>0)")" = "True" ]; then
    SELF_PCT=$(python3 -c "print(f'{float(\"$SELFING_RATE\")*100:g}')")
    SELF_TAG="_self${SELF_PCT}"; SELF_FLAG="--selfing-rate $SELFING_RATE"
fi

WORK=$CTRL/sims/cov${COVERAGE}_n${N_INDIV}_g${N_GEN}_s${SEED}${SELF_TAG}_hotspots_p80_chr1
mkdir -p $WORK $CTRL/logs

CACTUS_DIR=$CTRL/fastas_80
FOUNDERS_META=$CTRL/data/kmer_pa_p80/kmer_pa_Chr1.meta.npz
CN_KMER_PREFIX=$CTRL/data/kmer_pa_p80/kmer_pa
CN_VAR=$CTRL/data/var_pa_p80.var_pa.npz
CN_VAR_META=$CTRL/data/var_pa_p80.meta.npz
for f in "$CACTUS_DIR" "$FOUNDERS_META" "$CN_VAR" "$CN_VAR_META"; do
    [ -e "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

# N_INDIV > 80 is allowed: gen-0 founder draw is multinomial-with-replacement,
# so each founder can be assigned to multiple individuals (their lineages diverge
# during recombination). With n=200 from 80 founders, every founder is represented
# ~2.5× in expectation — the "perfect mix" regime analog of n200_g1 on the 231-panel.
if [ "$N_INDIV" -gt 80 ]; then
    echo "[$(date)] N_INDIV=$N_INDIV > 80 founders: with-replacement multinomial sampling"
fi

# -----------------------------------------------------------------------------
# STAGE 1: mosaic FASTAs
# -----------------------------------------------------------------------------
# SEEDMIX mimicry: always use no-replace / balanced allocation at gen-0 so the
# chrom-average truth_h per founder is uniform 1/min(n,F). Applies to g0 (truth
# is uniform per-founder) and g>=1 (recomb permutes ancestry within individuals
# but the chrom-average per-founder count is preserved by the gen-0 draw).
GEN0_FLAG="--gen0-no-replace"

echo "[$(date)] STAGE 1: make mosaic FASTAs (--chroms $CHROMS) $GEN0_FLAG ${SELF_FLAG:-(outcross)}"
$PYTHON /global/scratch/users/tbellg/kmate/sims/scripts/make_recomb_mosaics.py \
    --n-indiv $N_INDIV \
    --n-generations $N_GEN \
    --seed $SEED \
    --cactus-dir $CACTUS_DIR \
    --founders-meta $FOUNDERS_META \
    --out-dir $WORK \
    --chroms "$CHROMS" \
    $GEN0_FLAG $SELF_FLAG

# -----------------------------------------------------------------------------
# STAGE 2: VISOR SHORtS at cov10x, Chr1 only
# -----------------------------------------------------------------------------
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

    # Snap `o` to 0.01%, then iteratively NUDGE the absorber `d` until
    # VISOR's exact sum order yields 100.0 in fp. Per memory
    # `feedback_visor_clonefraction_fp`: VISOR uses Python's `sum([d, o, o, ...])`
    # which is left-to-right serial addition. fp non-associativity means the
    # residual depends on the addition order. Compute `d` via fixed-point
    # iteration on VISOR's actual sum order — converges in 1-3 iterations.
    # Use repr() to preserve full fp precision through bash and VISOR's argparse.
    read OTHER_FRAC FIRST_FRAC <<<"$(python3 -c "
N = $N_INDIV
o = round(100.0 / N * 100) / 100.0          # snap o to 0.01%
d = 100.0 - (N - 1) * o                      # initial estimate
for _trial in range(50):
    s = sum([d] + [o] * (N - 1))             # VISOR's exact sum semantics
    if s == 100.0:
        break
    d += 100.0 - s                           # nudge by residual
else:
    raise AssertionError(f'fp fixed-point did not converge after 50 iters; last sum={s!r}')
print(repr(o), repr(d))
")"
    CLONE_DIRS=()
    FRACS=()
    HAP_DIRS=( $WORK/haps/s_ind* )
    for k in "${!HAP_DIRS[@]}"; do
        CLONE_DIRS+=("${HAP_DIRS[$k]}")
        if [ "$k" = "0" ]; then
            FRACS+=("$FIRST_FRAC")
        else
            FRACS+=("$OTHER_FRAC")
        fi
    done
    echo "  ${#CLONE_DIRS[@]} mosaic clones, others=${OTHER_FRAC}% each, ind001=${FIRST_FRAC}% (absorbs rounding residual)"
    SUM_CHK=$(echo "${FRACS[@]}" | awk '{for(i=1;i<=NF;i++) s+=$i; printf "%.6f", s}')
    echo "  sum check: $SUM_CHK% (must be exactly 100)"

    source /global/home/users/tbellg/miniforge3/etc/profile.d/conda.sh && conda activate kmate
    rm -rf $READS_DIR; mkdir -p $READS_DIR
    VISOR SHORtS \
        -g $REF \
        -s "${CLONE_DIRS[@]}" \
        -b $REGION_BED \
        -o $READS_DIR \
        --coverage $COVERAGE \
        --clonefraction "${FRACS[@]}" \
        --error 0.001 \
        --length 150 \
        --fastq \
        --threads 8 || true
    if [ ! -s ${READS_DIR}/r1.fq ] || [ ! -s ${READS_DIR}/r2.fq ]; then
        echo "ERROR: VISOR SHORtS failed" >&2; exit 1
    fi
fi
echo "[$(date)] reads: $(du -h ${READS_DIR}/r1.fq ${READS_DIR}/r2.fq | tail -2)"

# -----------------------------------------------------------------------------
# STAGE 3: truth from ancestry tracks (var_pa_p80)
# -----------------------------------------------------------------------------
echo
echo "[$(date)] STAGE 3: per-record truth from ancestry tracks (var_pa_p80)"
$PYTHON /global/scratch/users/tbellg/kmate/sims/scripts/compute_recomb_truth.py \
    --ancestry $WORK/ancestry.tsv \
    --weights $WORK/pool_weights.tsv \
    --var-pa $CN_VAR \
    --var-meta $CN_VAR_META \
    --out $WORK/recomb_truth.tsv.gz

echo
echo "[$(date)] DONE -- $WORK"
ls -lh $WORK/reads/r1.fq $WORK/reads/r2.fq $WORK/recomb_truth.tsv.gz
