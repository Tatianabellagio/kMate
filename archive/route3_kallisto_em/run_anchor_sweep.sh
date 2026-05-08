#!/bin/bash
# Anchor-weight sweep for kallisto-EM. Reads from a previously saved EC cache,
# runs per-window EM with several global-anchor weights, evaluates each.
#
# Usage:
#   bash run_anchor_sweep.sh <ec_cache> <sim_name> <out_dir>
#
# Example:
#   bash run_anchor_sweep.sh /tmp/kallisto_em_full_ec_cache.npz \\
#                            cov50_n50_g3_s42_hotspots_p231_chr1 \\
#                            /tmp/anchor_sweep
set -euo pipefail
EC_CACHE=${1:?Usage: $0 <ec_cache> <sim_name> <out_dir>}
SIM=${2:?Usage: $0 <ec_cache> <sim_name> <out_dir>}
OUT_DIR=${3:?Usage: $0 <ec_cache> <sim_name> <out_dir>}
mkdir -p "$OUT_DIR"

PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python
SRC=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/src
SCRIPTS=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/scripts

# Anchor weights to test. 0 = no anchor (matches the original full-cov run).
ANCHORS=(0.0 0.05 0.1 0.3 1.0)

# Reads arg is required by argparse but unused when --load-ec-cache is set.
# Pass any FASTQ file just to satisfy the parser.
READS_R1=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/sims/visor_freqk/pool_sweep_82_recomb/${SIM}/reads/r1.fq

for L in "${ANCHORS[@]}"; do
    OUT_TSV="$OUT_DIR/kallisto_em_anchor${L}.tsv"
    OUT_H="$OUT_DIR/kallisto_em_anchor${L}_h.npz"
    LOG="$OUT_DIR/anchor${L}.log"
    if [[ -s $OUT_TSV ]]; then
        echo "[skip] $OUT_TSV exists"
        continue
    fi
    echo "=== anchor=$L → $OUT_TSV ==="
    $PYTHON -u $SRC/per_sample_kallisto_em.py \
        --cn-kmer-prefix /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_full_231_v2/cn \
        --cn-var /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_var_231_v2.cn_var.npz \
        --cn-var-meta /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_var_231_v2.meta.npz \
        --reads "$READS_R1" \
        --sample anchor${L} \
        --out "$OUT_TSV" \
        --chroms Chr1 --window-bp 10000 \
        --load-ec-cache "$EC_CACHE" \
        --save-h-blocks "$OUT_H" \
        --global-anchor-weight $L > "$LOG" 2>&1
    tail -5 "$LOG"
    echo
done

echo "=== sweep done. Evaluating ==="
for L in "${ANCHORS[@]}"; do
    OUT_TSV="$OUT_DIR/kallisto_em_anchor${L}.tsv"
    [[ -s $OUT_TSV ]] || continue
    echo "--- anchor=$L ---"
    $PYTHON $SCRIPTS/eval_kallisto_em.py \
        --sim "$SIM" --kal-tsv "$OUT_TSV" \
        --label "anchor${L}" 2>&1 | grep -E "^[ -]|R²" | head -20
done
