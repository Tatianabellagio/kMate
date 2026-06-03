#!/bin/bash
#SBATCH --job-name=jf_pair
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --array=1-82%50
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --requeue
#SBATCH --output=logs/pair_%A_%a.out
#SBATCH --error=logs/pair_%A_%a.err

# Each task: take anchor sample i (by array index), compute Jaccard vs samples i+1..82
# This produces upper-triangle rows of the 82x82 matrix
mkdir -p logs
set -euo pipefail
ROOT=/global/scratch/users/tbellg/hapfire_sv/jf_chr1
JF=/global/home/users/tbellg/miniforge3/envs/pangenie/bin/jellyfish

# Map array task → anchor sample
ANCHOR_LINE=${SLURM_ARRAY_TASK_ID}
mapfile -t SAMPLES < $ROOT/scripts/asm_ids.txt
NS=${#SAMPLES[@]}
if [ "$ANCHOR_LINE" -gt "$NS" ]; then
  echo "task $ANCHOR_LINE > $NS samples, exiting"; exit 0
fi
A=${SAMPLES[$((ANCHOR_LINE - 1))]}
OUT=$ROOT/pairs/${ANCHOR_LINE}_${A}.tsv

# Idempotent: if all pairs already written, skip
if [ -s "$OUT" ]; then
  N_EXPECT=$((NS - ANCHOR_LINE))
  N_HAVE=$(tail -n +2 "$OUT" | wc -l)
  if [ "$N_HAVE" -eq "$N_EXPECT" ]; then
    echo "[$(date +%T)] $A: $N_HAVE/$N_EXPECT pairs already done"
    exit 0
  fi
fi

# Header
echo -e "A\tB\tnA\tnB\tnU\tnI\tjaccard" > "$OUT"

# Cache |A|
[ -s "$ROOT/jf/$A.jf" ] || { echo "ERROR: $ROOT/jf/$A.jf missing"; exit 1; }
nA=$($JF stats $ROOT/jf/$A.jf | awk '/Distinct/{print $2}')

# Each pair: merge → |A∪B|, compute |A∩B| and Jaccard
TMP=$(mktemp -d -p /tmp jfpair_${A}_XXXX)
trap "rm -rf $TMP" EXIT

for ((j=ANCHOR_LINE; j<NS; j++)); do
  B=${SAMPLES[$j]}
  [ -s "$ROOT/jf/$B.jf" ] || { echo "WARN: missing $B.jf, skipping"; continue; }
  nB=$($JF stats $ROOT/jf/$B.jf | awk '/Distinct/{print $2}')
  MJF=$TMP/m.jf
  $JF merge -o $MJF $ROOT/jf/$A.jf $ROOT/jf/$B.jf 2>/dev/null
  nU=$($JF stats $MJF | awk '/Distinct/{print $2}')
  nI=$((nA + nB - nU))
  jac=$(awk "BEGIN{printf \"%.6f\", $nI/$nU}")
  echo -e "$A\t$B\t$nA\t$nB\t$nU\t$nI\t$jac" >> "$OUT"
  rm -f $MJF
done

echo "[$(date +%T)] $A done, $(tail -n +2 $OUT | wc -l) pairs written"
