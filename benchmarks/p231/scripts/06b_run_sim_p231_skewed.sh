#!/bin/bash
#SBATCH --job-name=p231_b_sim_skew
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH --output=logs/06b_sim_skew_%j.out
#SBATCH --error=logs/06b_sim_skew_%j.err

# =============================================================================
# benchmarks/p231 Phase B' (skewed) -- cov10 sim on the 231 panel with ONE DOMINANT
# individual (ind001 gets DOMINANT_FRAC% of pool reads). Mirrors benchmarks/p80
# 06b with the SAME two p231 changes: (1) RANDOM crossovers (no LD-block flag),
# (2) truth on BOTH arch3 var_pas (atomized + raw).
#
# Usage: sbatch 06b_run_sim_p231_skewed.sh N_INDIV N_GEN [SEED=42] [DOMINANT_FRAC=50.0]
#   for n50_g3_dom500 use: 50 3 42 50.0  (dir tag dom500)
# =============================================================================
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"; export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

N_INDIV=${1:?Usage: N_INDIV N_GEN [SEED=42] [DOMINANT_FRAC=50.0]}
N_GEN=${2:?Usage: N_INDIV N_GEN [SEED=42] [DOMINANT_FRAC=50.0]}
SEED=${3:-42}
DOMINANT_FRAC=${4:-50.0}
WINNER=${5:-recomb}   # recomb (default: don't force the winner's type) | pure (force non-recombinant winner)
SELFING_RATE=${6:-0.0}   # background selfing rate (0.97 => winner most-probably non-recombinant, no forcing)
COVERAGE=10; CHROMS="Chr1"

ROOT=/global/scratch/users/tbellg/kmate
CTRL=$ROOT/benchmarks/p231
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
PYTHON=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
SCRIPTS=$CTRL/scripts

DOM_TAG="dom$(echo $DOMINANT_FRAC | tr -d '.')"
PURE_FLAG=""
if [ "$WINNER" = "pure" ]; then DOM_TAG="${DOM_TAG}nr"; PURE_FLAG="--dominant-pure-founder"; fi
SELF_TAG=""; SELF_FLAG=""
if [ "$(python3 -c "print(float('$SELFING_RATE')>0)")" = "True" ]; then
    SELF_PCT=$(python3 -c "print(f'{float(\"$SELFING_RATE\")*100:g}')")
    SELF_TAG="_self${SELF_PCT}"; SELF_FLAG="--selfing-rate $SELFING_RATE"
fi
WORK=$CTRL/sims/cov${COVERAGE}_n${N_INDIV}_g${N_GEN}_s${SEED}${SELF_TAG}_hotspots_${DOM_TAG}_p231_chr1
mkdir -p $WORK $CTRL/logs

CACTUS_DIR=$CTRL/fastas_231
FOUNDERS_META=$ROOT/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1.meta.npz
CN_VAR_ATOM=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1_atomized.var_pa.npz
CN_VAR_ATOM_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1_atomized.meta.npz
CN_VAR_RAW=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
CN_VAR_RAW_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz
TRUTH=$ROOT/sims/scripts/compute_recomb_truth.py

for f in "$CACTUS_DIR" "$FOUNDERS_META" "$CN_VAR_ATOM" "$CN_VAR_RAW" "$TRUTH"; do
    [ -e "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done
echo "[$(date)] DOMINANT_FRAC=${DOMINANT_FRAC}% (ind001)"

# STAGE 1: mosaics -- RANDOM crossovers (NO --crossovers-from-ld-blocks)
echo "[$(date)] STAGE 1: mosaics (RANDOM crossovers @ 4 cM/Mb)"
$PYTHON $ROOT/sims/scripts/make_recomb_mosaics.py \
    --n-indiv $N_INDIV --n-generations $N_GEN --seed $SEED \
    --cactus-dir $CACTUS_DIR --founders-meta $FOUNDERS_META \
    --out-dir $WORK --chroms "$CHROMS" $PURE_FLAG $SELF_FLAG

# STAGE 2: VISOR with SKEWED fractions
echo; echo "[$(date)] STAGE 2: VISOR SHORtS cov${COVERAGE}x (ind001=${DOMINANT_FRAC}%)"
READS_DIR=$WORK/reads
if [ ! -s ${READS_DIR}/r1.fq ] || [ ! -s ${READS_DIR}/r2.fq ]; then
    REGION_BED=$WORK/region.bed
    declare -A CHROM_LEN=( [Chr1]=30427671 [Chr2]=19698289 [Chr3]=23459830 [Chr4]=18585056 [Chr5]=26975502 )
    > $REGION_BED
    for c in $CHROMS; do printf "%s\t1\t%d\t100.0\t100.0\n" "$c" "${CHROM_LEN[$c]}" >> $REGION_BED; done

    read OTHER_FRAC DOM_ACTUAL <<<"$(python3 -c "
o = round((100.0 - $DOMINANT_FRAC) / ($N_INDIV - 1) * 10) / 10.0
d = round(100.0 - ($N_INDIV - 1) * o, 1)
serial = d
for _ in range($N_INDIV - 1):
    serial += o
assert serial == 100.0, f'serial={serial!r} (snap-to-0.1 failed)'
print(repr(o), repr(d))
")"
    CLONE_DIRS=(); FRACS=(); HAP_DIRS=( $WORK/haps/s_ind* )
    for k in "${!HAP_DIRS[@]}"; do
        CLONE_DIRS+=("${HAP_DIRS[$k]}")
        if [ "$k" = "0" ]; then FRACS+=("$DOM_ACTUAL"); else FRACS+=("$OTHER_FRAC"); fi
    done
    echo "  ind001=${DOM_ACTUAL}% (target ${DOMINANT_FRAC}%), others=${OTHER_FRAC}% each"
    SUM_CHK=$(echo "${FRACS[@]}" | awk '{for(i=1;i<=NF;i++) s+=$i; printf "%.4f", s}')
    echo "  sum check: $SUM_CHK%"
    { echo -e "ind_id\tclone_dir\tvisor_pct"
      for k in "${!CLONE_DIRS[@]}"; do
        printf "%s\t%s\t%s\n" "$(basename "${CLONE_DIRS[$k]}")" "${CLONE_DIRS[$k]}" "${FRACS[$k]}"
      done; } > $WORK/visor_pool_fractions.tsv

    source /global/home/users/tbellg/miniforge3/etc/profile.d/conda.sh && conda activate kmate
    rm -rf $READS_DIR; mkdir -p $READS_DIR
    VISOR SHORtS -g $REF -s "${CLONE_DIRS[@]}" -b $REGION_BED -o $READS_DIR \
        --coverage $COVERAGE --clonefraction "${FRACS[@]}" \
        --error 0.001 --length 150 --fastq --threads 8 || true
    [ -s ${READS_DIR}/r1.fq ] && [ -s ${READS_DIR}/r2.fq ] || { echo "ERROR: VISOR failed" >&2; exit 1; }
    conda deactivate || true
fi
echo "[$(date)] reads: $(du -h ${READS_DIR}/r1.fq ${READS_DIR}/r2.fq | tail -2)"

# STAGE 3: truth on BOTH var_pas, using the SKEWED visor_pool_fractions.tsv
echo; echo "[$(date)] STAGE 3a: truth ATOMIZED (skewed weights)"
$PYTHON $TRUTH --ancestry $WORK/ancestry.tsv --weights $WORK/visor_pool_fractions.tsv \
    --var-pa $CN_VAR_ATOM --var-meta $CN_VAR_ATOM_META --out $WORK/recomb_truth_atomized.tsv.gz
echo; echo "[$(date)] STAGE 3b: truth RAW (skewed weights)"
$PYTHON $TRUTH --ancestry $WORK/ancestry.tsv --weights $WORK/visor_pool_fractions.tsv \
    --var-pa $CN_VAR_RAW --var-meta $CN_VAR_RAW_META --out $WORK/recomb_truth_raw.tsv.gz

echo; echo "[$(date)] DONE -- $WORK"
ls -lh $WORK/reads/r1.fq $WORK/recomb_truth_atomized.tsv.gz $WORK/recomb_truth_raw.tsv.gz
