#!/bin/bash
#SBATCH --job-name=qc_prep
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --output=logs/qc_%A_%a.out
#SBATCH --error=logs/qc_%A_%a.err

# =============================================================================
# qc_one.sh — per-ecotype sanity check on a preprocessed dedup'd fastq pair.
#
# Reads from preprocess_qc/output/qc_index_<panel>.tsv (built once by
# build_qc_index.py), one row per array task. Per row:
#   - gzip -t each file (CRC integrity)
#   - count reads via `zcat | wc -l / 4`
#   - sum bases via `zcat | awk '{if(NR%4==2) b+=length($0)} END{print b}'`
#   - PE: assert R1_reads == R2_reads
#
# Emits a single TSV row per task to preprocess_qc/output/qc_results/.
# =============================================================================
set -eo pipefail

BASE=/global/scratch/users/tbellg/hapfire_sv/preprocess_qc
INDEX=$BASE/output/qc_index_${PANEL:?must set PANEL=main|loo}.tsv
OUT_DIR=$BASE/output/qc_results
mkdir -p $OUT_DIR $BASE/logs

IDX=${SLURM_ARRAY_TASK_ID:-1}
LINE=$(sed -n "$((IDX+1))p" $INDEX)
[ -n "$LINE" ] || { echo "ERROR: no row at $IDX" >&2; exit 1; }

# index columns: panel ecotype layout file_r1 file_r2
PANEL_R=$(echo "$LINE" | cut -f1)
ECOTYPE=$(echo "$LINE" | cut -f2)
LAYOUT=$(echo "$LINE"  | cut -f3)
F1=$(echo "$LINE"      | cut -f4)
F2=$(echo "$LINE"      | cut -f5)

OUT=$OUT_DIR/${PANEL_R}_${ECOTYPE}.tsv

# function: produce "n_records<TAB>n_bases<TAB>min_len<TAB>max_len<TAB>mean_len" for a fq.gz
stats() {
    local f=$1
    zcat "$f" | awk '
        NR % 4 == 2 {
            n++; bases += length($0); l = length($0)
            if (n == 1) { min = max = l }
            if (l < min) min = l
            if (l > max) max = l
        }
        END {
            mean = (n>0) ? bases/n : 0
            printf "%d\t%d\t%d\t%d\t%.1f", n, bases, min, max, mean
        }'
}

# function: gzip integrity test; "OK" or "BAD"
gz_test() {
    if gzip -t "$1" 2>/dev/null; then echo OK; else echo BAD; fi
}

echo "[$(date)] QC ecotype=$ECOTYPE panel=$PANEL_R layout=$LAYOUT"
GZ1=$(gz_test "$F1")
S1=$(stats "$F1")
N1=$(echo "$S1" | cut -f1)
B1=$(echo "$S1" | cut -f2)
LMIN1=$(echo "$S1" | cut -f3)
LMAX1=$(echo "$S1" | cut -f4)
LMEAN1=$(echo "$S1" | cut -f5)

if [ "$LAYOUT" = "PE" ]; then
    GZ2=$(gz_test "$F2")
    S2=$(stats "$F2")
    N2=$(echo "$S2" | cut -f1)
    B2=$(echo "$S2" | cut -f2)
    LMIN2=$(echo "$S2" | cut -f3)
    LMAX2=$(echo "$S2" | cut -f4)
    LMEAN2=$(echo "$S2" | cut -f5)
    PAIR_OK="$([ "$N1" = "$N2" ] && echo yes || echo no)"
else
    GZ2="NA"; N2=0; B2=0; LMIN2=NA; LMAX2=NA; LMEAN2=NA; PAIR_OK="NA"
fi

# coverage estimate against TAIR10 chr-only ~119.7Mb (Chr1+...+Chr5)
TOTAL_B=$((B1 + B2))
COV=$(awk -v b=$TOTAL_B 'BEGIN{printf "%.2f", b/119700000}')

{
    echo "panel	ecotype	layout	gz_r1	gz_r2	n_r1	n_r2	pair_ok	bases_total	cov_est	lmin_r1	lmax_r1	lmean_r1	lmin_r2	lmax_r2	lmean_r2"
    echo "$PANEL_R	$ECOTYPE	$LAYOUT	$GZ1	$GZ2	$N1	$N2	$PAIR_OK	$TOTAL_B	$COV	$LMIN1	$LMAX1	$LMEAN1	$LMIN2	$LMAX2	$LMEAN2"
} > $OUT

cat $OUT
echo "[$(date)] DONE $ECOTYPE"
