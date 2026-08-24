#!/usr/bin/env bash
# Push seedmix (gen0) + gen1/2/3 allele-frequency matrices (snp/smallindel/sv) to
# the MOI-LAB shared Drive so the advisor can run analysis on them directly.
#
# data/class_matrices/      <- phase1_replication/results/class_matrices (gen0/1/2/3 only):
#               MAF-filtered, LD-block-tagged, true snp/sv/smallindel split.
#               Rebuilt 2026-07-27 on the post-Kf_w-fix pool_matrices (the prior
#               gen1/gen3 files there were stale, pre-fix; gen2 didn't exist).
# data-raw/gen_matrices/, data-raw/seedmix/  <- gen_matrices (gen1/2/3, unfiltered,
#               per-sample; snp/nonsnp/smallindel -- nonsnp = indel+SV combined,
#               no clean sv-only file) + the 8 raw SEEDMIX rep TSVs (founding,
#               genome-wide, unsplit).
#
# Run (survives logout):
#   nohup bash analysis/grenenet_selection/push_af_matrices_to_drive.sh > /tmp/af_push.log 2>&1 &
#   tail -f /tmp/af_push.log
# PREREQ: valid rclone token. Verify: rclone lsd gdrive:PROJECTS/grenenet/GrENE-net_PHASE1SV/
set -euo pipefail

PROJ=/global/scratch/projects/fc_moilab/tbellg/kmate
GEA="$PROJ/analysis/grenenet_selection"
DEST_BASE="gdrive:PROJECTS/grenenet/GrENE-net_PHASE1SV"
RCLONE="rclone"
FLAGS=(--stats 30s --stats-one-line -v --transfers 8 --checkers 16
       --timeout 120s --contimeout 30s --low-level-retries 20 --retries 10 --retries-sleep 10s)

echo "### started: $(date)"

echo "### [1/3] processed class_matrices (gen0/1/2/3, snp/sv/smallindel) -> data/class_matrices/"
for CLS in snp sv smallindel; do
  for G in 0 1 2 3; do
    $RCLONE copy "$GEA/r2_gea_nonsnp/phase1_replication/results/class_matrices" \
      "$DEST_BASE/data/class_matrices" \
      --include "${CLS}_gen${G}_af.npy" --include "${CLS}_gen${G}.records.csv" "${FLAGS[@]}"
  done
done
$RCLONE copy "$GEA/r2_gea_nonsnp/phase1_replication/results/class_matrices" \
  "$DEST_BASE/data/class_matrices" \
  --include "gen{0,1,2,3}.pools.csv" "${FLAGS[@]}"

echo "### [2/3] raw gen_matrices (gen1/2/3, snp/nonsnp/smallindel) -> data-raw/gen_matrices/"
$RCLONE copy "$GEA/common/results/gen_matrices" "$DEST_BASE/data-raw/gen_matrices" \
  --include "gen{1,2,3}_{snp,nonsnp,smallindel}_af.npy" \
  --include "gen{1,2,3}.rowmeta.csv" "${FLAGS[@]}"

echo "### [3/3] raw SEEDMIX founding reps (genome-wide, unsplit) -> data-raw/seedmix/"
$RCLONE copy "$GEA/common/rerun_kfw_hb/seedmix" "$DEST_BASE/data-raw/seedmix" \
  --include "SEEDMIX_S[1-8].tsv" "${FLAGS[@]}"

echo "### finished: $(date)"
echo "### verify:  rclone ls \"$DEST_BASE/data\" ; rclone ls \"$DEST_BASE/data-raw\""
