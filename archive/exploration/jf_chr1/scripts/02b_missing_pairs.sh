#!/bin/bash
#SBATCH --job-name=jf_miss
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --array=1-120%80
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=00:10:00
#SBATCH --requeue
#SBATCH --output=logs/miss_%A_%a.out
#SBATCH --error=logs/miss_%A_%a.err

mkdir -p logs
set -euo pipefail
ROOT=/global/scratch/users/tbellg/hapfire_sv/jf_chr1
JF=/global/home/users/tbellg/miniforge3/envs/pangenie/bin/jellyfish

# Pull this task's pair
PAIR_LINE=${SLURM_ARRAY_TASK_ID}
LINE=$(sed -n "${PAIR_LINE}p" $ROOT/scripts/missing_pairs.txt)
A=$(echo "$LINE" | cut -f1)
B=$(echo "$LINE" | cut -f2)
OUT_DIR=$ROOT/pairs_missing
mkdir -p $OUT_DIR
OUT=$OUT_DIR/${A}_${B}.tsv
[ -s "$OUT" ] && { echo "already done"; exit 0; }

TMP=$(mktemp -d -p /tmp jfm_${A}_${B}_XXXX)
trap "rm -rf $TMP" EXIT

nA=$($JF stats $ROOT/jf/$A.jf | awk '/Distinct/{print $2}')
nB=$($JF stats $ROOT/jf/$B.jf | awk '/Distinct/{print $2}')
$JF merge -o $TMP/m.jf $ROOT/jf/$A.jf $ROOT/jf/$B.jf 2>/dev/null
nU=$($JF stats $TMP/m.jf | awk '/Distinct/{print $2}')
nI=$((nA + nB - nU))
jac=$(awk "BEGIN{printf \"%.6f\", $nI/$nU}")
echo -e "A\tB\tnA\tnB\tnU\tnI\tjaccard" > "$OUT"
echo -e "$A\t$B\t$nA\t$nB\t$nU\t$nI\t$jac" >> "$OUT"
echo "[$(date +%T)] $A $B  J=$jac"
