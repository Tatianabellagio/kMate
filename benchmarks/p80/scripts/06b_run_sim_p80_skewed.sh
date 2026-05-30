#!/bin/bash
#SBATCH --job-name=p80_b_sim_skew
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/06b_sim_skew_%j.out
#SBATCH --error=logs/06b_sim_skew_%j.err

# =============================================================================
# Phase B' (skewed) -- recombinant cov10 sim on p80 with ONE DOMINANT individual.
#
# Same mosaic-generation as 06_run_sim_p80.sh, but VISOR pool fractions are
# SKEWED: one mosaic (ind001) gets DOMINANT_FRAC% of pool reads; the other
# N_INDIV-1 mosaics split the remainder equally.
#
# Models "selection-like" stress: a single ancestry pattern dominates the
# chrom-wide signal, breaking ergodicity → per-position founder mix becomes
# spatially non-stationary. Stress-tests global vs window methods.
#
# Usage:
#   sbatch 06b_run_sim_p80_skewed.sh N_INDIV N_GEN [SEED=42] [DOMINANT_FRAC=50.0]
# =============================================================================
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

N_INDIV=${1:?Usage: sbatch 06b_run_sim_p80_skewed.sh N_INDIV N_GEN [SEED=42] [DOMINANT_FRAC=50.0]}
N_GEN=${2:?Usage: sbatch 06b_run_sim_p80_skewed.sh N_INDIV N_GEN [SEED=42] [DOMINANT_FRAC=50.0]}
SEED=${3:-42}
DOMINANT_FRAC=${4:-50.0}
COVERAGE=10
CHROMS="Chr1"

CTRL=/global/scratch/users/tbellg/kmate/benchmarks/p80
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
SCRIPTS=$CTRL/scripts

DOM_TAG="dom$(echo $DOMINANT_FRAC | tr -d '.')"
WORK=$CTRL/sims/cov${COVERAGE}_n${N_INDIV}_g${N_GEN}_s${SEED}_hotspots_${DOM_TAG}_p80_chr1
mkdir -p $WORK $CTRL/logs

CACTUS_DIR=$CTRL/fastas_80
FOUNDERS_META=$CTRL/data/kmer_pa_p80/kmer_pa_Chr1.meta.npz
CN_VAR=$CTRL/data/var_pa_p80.var_pa.npz
CN_VAR_META=$CTRL/data/var_pa_p80.meta.npz
for f in "$CACTUS_DIR" "$FOUNDERS_META" "$CN_VAR" "$CN_VAR_META"; do
    [ -e "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

if [ "$N_INDIV" -gt 80 ]; then
    echo "[$(date)] N_INDIV=$N_INDIV > 80 founders: with-replacement multinomial sampling"
fi
echo "[$(date)] DOMINANT_FRAC=${DOMINANT_FRAC}% (ind001), others split $(echo "scale=2; (100 - $DOMINANT_FRAC) / ($N_INDIV - 1)" | bc)% each"

# -----------------------------------------------------------------------------
# STAGE 1: mosaic FASTAs (identical to non-skewed)
# -----------------------------------------------------------------------------
echo "[$(date)] STAGE 1: make mosaic FASTAs (--chroms $CHROMS)"
$PYTHON /global/scratch/users/tbellg/kmate/sims/scripts/make_recomb_mosaics.py \
    --n-indiv $N_INDIV \
    --n-generations $N_GEN \
    --seed $SEED \
    --cactus-dir $CACTUS_DIR \
    --founders-meta $FOUNDERS_META \
    --out-dir $WORK \
    --chroms "$CHROMS"

# -----------------------------------------------------------------------------
# STAGE 2: VISOR SHORtS with SKEWED fractions
# -----------------------------------------------------------------------------
echo
echo "[$(date)] STAGE 2: VISOR SHORtS at cov ${COVERAGE}x (skewed: ind001=${DOMINANT_FRAC}%)"
READS_DIR=$WORK/reads
if [ ! -s ${READS_DIR}/r1.fq ] || [ ! -s ${READS_DIR}/r2.fq ]; then
    REGION_BED=$WORK/region.bed
    declare -A CHROM_LEN=( [Chr1]=30427671 [Chr2]=19698289 [Chr3]=23459830 [Chr4]=18585056 [Chr5]=26975502 )
    > $REGION_BED
    for c in $CHROMS; do
        printf "%s\t1\t%d\t100.0\t100.0\n" "$c" "${CHROM_LEN[$c]}" >> $REGION_BED
    done

    # Snap-to-0.1% pattern: VISOR requires strict-fp sum==100.0. The dominant
    # absorbs the rounding residual. See memory feedback_visor_clonefraction_fp.
    read OTHER_FRAC DOM_ACTUAL <<<"$(python3 -c "
o = round((100.0 - $DOMINANT_FRAC) / ($N_INDIV - 1) * 10) / 10.0
d = round(100.0 - ($N_INDIV - 1) * o, 1)
serial = d
for _ in range($N_INDIV - 1):
    serial += o
assert serial == 100.0, f'serial={serial!r} (snap-to-0.1 failed)'
print(repr(o), repr(d))
")"

    CLONE_DIRS=()
    FRACS=()
    HAP_DIRS=( $WORK/haps/s_ind* )
    for k in "${!HAP_DIRS[@]}"; do
        CLONE_DIRS+=("${HAP_DIRS[$k]}")
        if [ "$k" = "0" ]; then
            FRACS+=("$DOM_ACTUAL")
        else
            FRACS+=("$OTHER_FRAC")
        fi
    done
    echo "  ind001 absorbs rounding: ${DOM_ACTUAL}% (target was ${DOMINANT_FRAC}%)"
    echo "  ${#CLONE_DIRS[@]} clones; ind001=${DOM_ACTUAL}%, others=${OTHER_FRAC}% each"
    SUM_CHK=$(echo "${FRACS[@]}" | awk '{for(i=1;i<=NF;i++) s+=$i; printf "%.4f", s}')
    echo "  sum check: $SUM_CHK% (should be 100)"

    # Persist per-individual visor fractions so compute_recomb_truth uses them.
    {
        echo -e "ind_id\tclone_dir\tvisor_pct"
        for k in "${!CLONE_DIRS[@]}"; do
            ind_id=$(basename "${CLONE_DIRS[$k]}")
            printf "%s\t%s\t%s\n" "$ind_id" "${CLONE_DIRS[$k]}" "${FRACS[$k]}"
        done
    } > $WORK/visor_pool_fractions.tsv

    source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang
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
# STAGE 3: truth (uses visor_pool_fractions.tsv → skewed-mode in
# compute_recomb_truth.py automatically activates).
# -----------------------------------------------------------------------------
echo
echo "[$(date)] STAGE 3: per-record truth (skewed visor fractions)"
$PYTHON /global/scratch/users/tbellg/kmate/sims/scripts/compute_recomb_truth.py \
    --ancestry $WORK/ancestry.tsv \
    --weights  $WORK/visor_pool_fractions.tsv \
    --var-pa   $CN_VAR \
    --var-meta $CN_VAR_META \
    --out      $WORK/recomb_truth.tsv.gz

echo
echo "[$(date)] DONE -- $WORK"
ls -lh $WORK/reads/r1.fq $WORK/reads/r2.fq $WORK/recomb_truth.tsv.gz
