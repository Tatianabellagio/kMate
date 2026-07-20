#!/bin/bash
# Overnight self-healing monitor for the per_founder rerun cohort (robust v2).
# Tracks cohort jobs by NAME (kmate_grenenet), tolerant of transient squeue errors,
# requires 2 consecutive confirmed-empty polls before declaring drained, and chunks
# requeues to <=1000. Requeues any sample missing its final <SAMPLE>.tsv at 48G mem
# (transient jellyfish/dev-shm/OOM failures), up to MAXROUND rounds. Reports tally.
set -uo pipefail
cd /global/scratch/users/tbellg/kmate
BASE=analysis/grenenet_gea/rerun_kfw_hb
EVO_MAN=data/sample_manifest_usesample.tsv
SM_MAN=data/seedmix_manifest_arch3.tsv
RUNNER=grenenet/run_site_array_perchrom.sh
NAME=kmate_grenenet
COMMON="--mem=48G --cpus-per-task=8 --time=3:00:00 --partition=savio4_htc --account=co_moilab --qos=savio_lowprio"
MAXROUND=6
mkdir -p logs/reconcile

wait_drain () {   # block until NO job named $NAME is queued/running (robust to transient squeue errors)
  local empties=0 out rc
  while true; do
    out=$(squeue -u "$USER" -h -n "$NAME" -o "%i" 2>/dev/null); rc=$?
    if [ $rc -ne 0 ]; then sleep 120; continue; fi          # transient scheduler error -> retry, do NOT declare drained
    if [ -z "${out//[[:space:]]/}" ]; then
      empties=$((empties+1)); [ $empties -ge 2 ] && return 0 # 2 consecutive confirmed-empty polls
    else empties=0; fi
    sleep 180
  done
}

build_missing () {  # $1=orig manifest $2=out dir $3=missing-manifest-out ; echoes count
  head -1 "$1" > "$3"
  tail -n +2 "$1" | while IFS= read -r line; do
    sid=$(printf '%s' "$line" | cut -f1); [ -n "$sid" ] || continue
    [ -s "$2/${sid}.tsv" ] || printf '%s\n' "$line"
  done >> "$3"
  echo $(( $(wc -l < "$3") - 1 ))
}

requeue () {  # $1=missing-manifest $2=out dir $3=round $4=kind ; submits <=1000-row chunks
  local MAN=$1 DIR=$2 rd=$3 kind=$4 n off N th any=0
  n=$(( $(wc -l < "$MAN") - 1 )); [ "$n" -gt 0 ] || { echo 0; return; }
  off=0
  while [ $off -lt $n ]; do
    N=$(( n - off )); [ $N -gt 1000 ] && N=1000
    th=$([ $N -lt 120 ] && echo $N || echo 120)
    sbatch --parsable $COMMON --array=1-${N}%${th} \
      --output=logs/reconcile/rq_${kind}_r${rd}_o${off}_%A_%a.log \
      --export=ALL,MANIFEST=$MAN,OUT_DIR=$DIR,BLOCK_MODE=global,KMER_WEIGHT=uniform,NORMALIZE=per_founder,OFFSET=$off \
      "$RUNNER" >/dev/null && any=1
    off=$(( off + 1000 ))
  done
  echo $any
}

for round in $(seq 1 $MAXROUND); do
  wait_drain
  echo "[round $round $(date)] cohort drained; scanning for missing"
  did=0
  for kind in sm evo; do
    if [ "$kind" = evo ]; then MAN=$EVO_MAN; DIR=$BASE/evolved; else MAN=$SM_MAN; DIR=$BASE/seedmix; fi
    MISS=logs/reconcile/missing_${kind}_r${round}.tsv
    n=$(build_missing "$MAN" "$DIR" "$MISS")
    echo "[round $round] $kind missing=$n"
    if [ "$n" -gt 0 ]; then a=$(requeue "$MISS" "$DIR" "$round" "$kind"); [ "$a" = 1 ] && did=1; fi
  done
  [ $did -eq 0 ] && { echo "[round $round $(date)] COHORT COMPLETE"; break; }
  sleep 60   # let the new arrays register before wait_drain re-polls
done

echo "=================== FINAL TALLY $(date) ==================="
echo "seed-mix final TSVs: $(ls $BASE/seedmix/*.tsv 2>/dev/null | grep -vcE '_Chr[0-9]') / 8"
echo "evolved  final TSVs: $(ls $BASE/evolved/*.tsv 2>/dev/null | grep -vcE '_Chr[0-9]') / 2168"
LM=logs/reconcile/_final_missing.txt
tail -n +2 "$EVO_MAN" | cut -f1 | while read s; do [ -s "$BASE/evolved/${s}.tsv" ] || echo "$s"; done > "$LM"
echo "evolved still-missing: $(grep -c . "$LM")  (list: $LM)"; head -20 "$LM"
echo "==========================================================="
