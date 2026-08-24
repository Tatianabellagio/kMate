#!/bin/bash
# ============================================================================
# Orchestrate the two-phase cohort runner (Phase A -> B -> C) in WAVES.
#
# ⚠️  OPTIONAL / NOT YET RUN AT SCALE. The validated production path is the
# single-phase grenenet/run_site_array_perchrom.sh (pilot: 175 samples,
# 2026-06-01). Two-phase trades simplicity for extra fan-out, but its per-sample
# DB lives on shared scratch (Phase A and Phase B can land on different nodes),
# so it is SUBJECT TO the "Vessel-B" shared-filesystem stall that the
# single-phase runner avoids by keeping the DB node-local. Prefer single-phase
# unless you specifically need the flat N*5 fan-out. See grenenet/README.md.
#
# Why waves: the cluster MaxArraySize is 1001, so no array can exceed 1000
# indices. A wave of W<=200 samples makes Phase B exactly W*5<=1000 tasks (one
# array) and caps live DB scratch at W*~5 GB (200 -> ~1 TB). Waves are chained
# A->B->C with SLURM dependencies; by default each wave's Phase A waits for the
# PREVIOUS wave's Phase C, so at most one wave's DBs sit on scratch at a time.
# (For more throughput at the cost of ~2x scratch, raise PIPELINE_DEPTH to 2.)
#
# Usage:
#   grenenet/submit_two_phase_cohort.sh MANIFEST OUT_DIR DB_DIR \
#       [WAVE_SIZE=200] [BLOCK_MODE=global] [THROTTLE=200]
#
# Example (full evolved cohort, --unit chrom — the production set):
#   grenenet/submit_two_phase_cohort.sh data/sample_manifest_usesample.tsv \
#       analysis/grenenet_gea/common/rerun_kfw_hb/evolved $SCRATCH/kmate_dbs 200 chrom 200
# ============================================================================
set -euo pipefail
cd /global/scratch/users/tbellg/kmate

MANIFEST=${1:?MANIFEST}; OUT_DIR=${2:?OUT_DIR}; DB_DIR=${3:?DB_DIR}
WAVE_SIZE=${4:-200}            # <=200 so Phase B array (WAVE_SIZE*5) <= 1000
BLOCK_MODE=${5:-global}
THROTTLE=${6:-200}             # %C concurrent-task cap on each array
PIPELINE_DEPTH=${PIPELINE_DEPTH:-1}   # waves allowed in flight (1 = min scratch)
NCHR=5

[ $((WAVE_SIZE * NCHR)) -le 1000 ] || { echo "WAVE_SIZE*5 must be <=1000 (MaxArraySize)"; exit 1; }
mkdir -p "$OUT_DIR" "$DB_DIR" logs

N=$(( $(wc -l < "$MANIFEST") - 1 ))   # minus header
echo "cohort: N=$N samples | wave=$WAVE_SIZE | mode=$BLOCK_MODE | throttle=$THROTTLE | pipeline_depth=$PIPELINE_DEPTH"

declare -a CJOBS=()   # Phase C job id per wave (for pipeline dependency)
w=0; off=0
while [ $off -lt $N ]; do
    w=$((w+1)); rem=$((N - off))
    W=$(( rem < WAVE_SIZE ? rem : WAVE_SIZE ))
    WB=$(( W * NCHR )); offB=$(( off * NCHR ))

    # Gate this wave's Phase A on the Phase C of the wave PIPELINE_DEPTH back, so
    # at most PIPELINE_DEPTH waves' DBs sit on scratch at once. CJOBS is 1-indexed
    # by wave number, so wave w gates on CJOBS[w - PIPELINE_DEPTH] (>= 1).
    depA=""
    gate_idx=$(( w - PIPELINE_DEPTH ))
    if [ $gate_idx -ge 1 ]; then depA="--dependency=afterok:${CJOBS[$gate_idx]}"; fi

    echo "=== wave $w: samples $((off+1))..$((off+W))  (OFFSET=$off) ${depA:+[gated on ${CJOBS[$gate_idx]}]} ==="

    jA=$(sbatch --parsable --array=1-${W}%${THROTTLE} $depA \
          --export=ALL,MANIFEST=$MANIFEST,DB_DIR=$DB_DIR,OFFSET=$off \
          grenenet/run_phaseA_build_dbs.sh)
    jB=$(sbatch --parsable --array=1-${WB}%${THROTTLE} --dependency=afterok:$jA \
          --export=ALL,MANIFEST=$MANIFEST,DB_DIR=$DB_DIR,OUT_DIR=$OUT_DIR,OFFSET=$offB,BLOCK_MODE=$BLOCK_MODE \
          grenenet/run_phaseB_em_per_chrom.sh)
    jC=$(sbatch --parsable --array=1-${W}%${THROTTLE} --dependency=afterok:$jB \
          --export=ALL,MANIFEST=$MANIFEST,DB_DIR=$DB_DIR,OUT_DIR=$OUT_DIR,OFFSET=$off \
          grenenet/run_phaseC_concat.sh)
    echo "  A=$jA  B=$jB  C=$jC"

    CJOBS[$w]=$jC
    off=$(( off + W ))
done
echo "submitted $w waves (A->B->C each). Watch: squeue -u $USER | grep kmate_phase"
