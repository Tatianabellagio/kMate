#!/bin/bash
#SBATCH --job-name=p231_b_sim
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH --output=logs/06_sim_%j.out
#SBATCH --error=logs/06_sim_%j.err

# =============================================================================
# control_p231 Phase B -- recombinant cov10 sim on the 231 panel (Chr1).
# Mirrors control_p80/06_run_sim_p80.sh with TWO deliberate changes:
#   (1) RANDOM crossover placement at the A. thaliana rate (4 cM/Mb): we DROP
#       --crossovers-from-ld-blocks so sample_crossovers() falls back to uniform
#       random positions. Forcing crossovers at hapFIRE BigLD boundaries was
#       circular for a benchmark (it biases toward the block-based method).
#   (2) per-record TRUTH built on BOTH arch3 cn_vars (atomized + raw), so each
#       projection arm (07c) joins its truth 100% on (chrom,pos,ref_len,alt_len).
#
# Usage: sbatch 06_run_sim_p231.sh N_INDIV N_GEN [SEED=42]
#   e.g. 50 0 ; 231 0 ; 50 1 ; 231 1 ; 50 3
# =============================================================================
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

N_INDIV=${1:?Usage: sbatch 06_run_sim_p231.sh N_INDIV N_GEN [SEED=42]}
N_GEN=${2:?Usage: sbatch 06_run_sim_p231.sh N_INDIV N_GEN [SEED=42]}
SEED=${3:-42}
COVERAGE=10
CHROMS="Chr1"

ROOT=/global/scratch/users/tbellg/kmate
CTRL=$ROOT/control_p231
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
SCRIPTS=$CTRL/scripts

WORK=$CTRL/sims/cov${COVERAGE}_n${N_INDIV}_g${N_GEN}_s${SEED}_hotspots_p231_chr1
mkdir -p $WORK $CTRL/logs

CACTUS_DIR=$CTRL/fastas_231
FOUNDERS_META=$ROOT/data/cn_full_231_v3qc_v3_filt2/cn_Chr1.meta.npz
# arch3 cn_vars (REUSED): atomized (SNP-level) + raw (SNP/indel/SV classes)
CN_VAR_ATOM=$ROOT/arch3/chr1/cn_var_231_arch3_chr1_atomized.cn_var.npz
CN_VAR_ATOM_META=$ROOT/arch3/chr1/cn_var_231_arch3_chr1_atomized.meta.npz
CN_VAR_RAW=$ROOT/arch3/chr1/cn_var_231_arch3_chr1.cn_var.npz
CN_VAR_RAW_META=$ROOT/arch3/chr1/cn_var_231_arch3_chr1.meta.npz
TRUTH=$ROOT/sims/visor_freqk/scripts/compute_recomb_truth.py

for f in "$CACTUS_DIR" "$FOUNDERS_META" "$CN_VAR_ATOM" "$CN_VAR_RAW" "$TRUTH"; do
    [ -e "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

# -----------------------------------------------------------------------------
# STAGE 1: mosaic FASTAs -- RANDOM crossovers (NO --crossovers-from-ld-blocks)
# -----------------------------------------------------------------------------
# SEEDMIX mimicry: always use no-replace / balanced allocation at gen-0 so the
# chrom-average truth_h per founder is uniform 1/min(n,F). Applies to g0 and
# g>=1 (recomb permutes ancestry within individuals but the per-founder draw
# count is preserved).
GEN0_FLAG="--gen0-no-replace"

echo "[$(date)] STAGE 1: mosaics (RANDOM crossovers @ 4 cM/Mb, --chroms $CHROMS) $GEN0_FLAG"
$PYTHON $SCRIPTS/make_recomb_mosaics_p231.py \
    --n-indiv $N_INDIV \
    --n-generations $N_GEN \
    --seed $SEED \
    --cactus-dir $CACTUS_DIR \
    --founders-meta $FOUNDERS_META \
    --out-dir $WORK \
    --chroms "$CHROMS" \
    $GEN0_FLAG

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

    # VISOR requires strict-fp clonefraction sum==100.0; snap o to 0.01% and
    # absorb residual in ind001 via fixed-point on VISOR's serial sum order
    # (memory feedback_visor_clonefraction_fp).
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
    SUM_CHK=$(echo "${FRACS[@]}" | awk '{for(i=1;i<=NF;i++) s+=$i; printf "%.6f", s}')
    echo "  sum check: $SUM_CHK% (must be exactly 100)"

    source "$(mamba info --base)/etc/profile.d/conda.sh" && conda activate pang
    rm -rf $READS_DIR; mkdir -p $READS_DIR
    VISOR SHORtS \
        -g $REF -s "${CLONE_DIRS[@]}" -b $REGION_BED -o $READS_DIR \
        --coverage $COVERAGE --clonefraction "${FRACS[@]}" \
        --error 0.001 --length 150 --fastq --threads 8 || true
    if [ ! -s ${READS_DIR}/r1.fq ] || [ ! -s ${READS_DIR}/r2.fq ]; then
        echo "ERROR: VISOR SHORtS failed" >&2; exit 1
    fi
    conda deactivate || true
fi
echo "[$(date)] reads: $(du -h ${READS_DIR}/r1.fq ${READS_DIR}/r2.fq | tail -2)"

# -----------------------------------------------------------------------------
# STAGE 3: per-record truth on BOTH cn_vars (atomized + raw)
# -----------------------------------------------------------------------------
echo
echo "[$(date)] STAGE 3a: truth on ATOMIZED cn_var (7.46M per-base records)"
$PYTHON $TRUTH \
    --ancestry $WORK/ancestry.tsv --weights $WORK/pool_weights.tsv \
    --cn-var $CN_VAR_ATOM --cn-var-meta $CN_VAR_ATOM_META \
    --out $WORK/recomb_truth_atomized.tsv.gz

echo
echo "[$(date)] STAGE 3b: truth on RAW arch3 cn_var (2.62M records, SNP/indel/SV)"
$PYTHON $TRUTH \
    --ancestry $WORK/ancestry.tsv --weights $WORK/pool_weights.tsv \
    --cn-var $CN_VAR_RAW --cn-var-meta $CN_VAR_RAW_META \
    --out $WORK/recomb_truth_raw.tsv.gz

echo
echo "[$(date)] DONE -- $WORK"
ls -lh $WORK/reads/r1.fq $WORK/reads/r2.fq $WORK/recomb_truth_atomized.tsv.gz $WORK/recomb_truth_raw.tsv.gz
