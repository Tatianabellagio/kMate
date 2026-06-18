#!/usr/bin/env bash
# Push ONLY the final GrENE-Net phase-1 GEA results to the MOI-LAB shared Drive,
# so the advisor can plot them himself for talks.
#
# Results-only by design:
#   - NO code           -> it lives on GitHub (Tatianabellagio/kMate)
#   - NO raw allele freq -> af_store (52 G) + gen/pool matrices (78 G) are huge AND
#                           rebuildable from the pipeline; they don't belong on Drive
#
# What goes up (~7 MB):
#   - phase-1 last-gen (gen9, bio1) WZA *per-block p-values*, primary deg7-cap2000
#     correction: 3 models (kendall / lfmm / binomial) x 3 classes (snp/smallindel/sv) = 9 CSVs
#   - the 3-model Manhattan figure (reference PNG)
#
# Run (survives logout):
#   nohup bash analysis/grenenet_gea/push_gea_to_drive.sh > /tmp/gea_push.log 2>&1 &
#   tail -f /tmp/gea_push.log
# PREREQ: valid rclone token.  Verify:  rclone lsd gdrive:PROJECTS/grenenet/
set -euo pipefail

PROJ=/global/scratch/users/tbellg/kmate
# Mirror the cluster folder names (grenenet_gea / phase1_replication / wza) under data/,
# but only rclone the selected result files.
DEST="gdrive:PROJECTS/grenenet/GrENE-net_PHASE1SV/data/grenenet_gea/phase1_replication"
RCLONE="rclone"
FLAGS=(--stats 10s --stats-one-line -v --transfers 8 --checkers 16
       --timeout 120s --contimeout 30s --low-level-retries 20 --retries 10 --retries-sleep 10s)

WZA="$PROJ/results/grenenet_gea/phase1_replication/wza"
P1="$PROJ/results/grenenet_gea/phase1_replication"

echo "### destination: $DEST"
echo "### started: $(date)"

# 1. Final WZA per-block p-values: gen9, bio1, deg7-cap2000, 3 models x 3 classes (9 files)
echo "### [1/2] phase-1 last-gen WZA per-block p-values -> wza/"
$RCLONE copy "$WZA" "$DEST/wza" \
    --include "wza_*_gen9_bio1_deg7cap2000.csv" "${FLAGS[@]}"

# 2. Reference figure (3-model Manhattan grid), PNG + PDF
echo "### [2/2] 3-model Manhattan figure (png + pdf) -> ./"
$RCLONE copy "$P1" "$DEST" \
    --include "manhattan_3models_deg7cap2000.png" \
    --include "manhattan_3models_deg7cap2000.pdf" "${FLAGS[@]}"

echo "### finished: $(date)"
echo "### verify:  rclone ls \"$DEST\""
