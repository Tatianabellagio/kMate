#!/usr/bin/env bash
# Push the GrENE-Net SV-GEA working set from the cluster to the MOI-LAB shared
# Drive, so it survives the week-long cluster outage. Server-side via rclone
# (gdrive: = MOI-LAB shared drive root). Run with nohup so it survives logout:
#
#   nohup bash analysis/grenenet_gea/push_gea_to_drive.sh > /tmp/gea_push.log 2>&1 &
#   tail -f /tmp/gea_push.log
#
# PREREQ: the rclone token must be valid. If expired, reconnect first (see README
# / the chat). Verify with:  rclone lsd gdrive:PROJECTS/grenenet/
set -euo pipefail

PROJ=/global/scratch/users/tbellg/kmate
DEST="gdrive:PROJECTS/grenenet/GrENE-net_PHASE1SV/kmate_gea_export"
RCLONE="rclone"
FLAGS=(--stats 20s --stats-one-line -v --transfers 8 --checkers 16
       --drive-chunk-size 128M --drive-acknowledge-abuse --fast-list
       # stall-resistance: drop hung sockets and retry instead of freezing
       --timeout 120s --contimeout 30s --expect-continue-timeout 30s
       --low-level-retries 20 --retries 10 --retries-sleep 10s)

echo "### destination: $DEST"
echo "### started: $(date)"

# 1. Code (the GEA pipeline itself)
echo "### [1/3] code -> code/"
$RCLONE copy "$PROJ/analysis/grenenet_gea" "$DEST/code" \
    --exclude "__pycache__/**" --exclude "*.pyc" "${FLAGS[@]}"

# 2. External inputs (scattered absolute paths the code reads; see MANIFEST.txt)
echo "### [2/3] external inputs -> external/"
declare -a EXT=(
  "/global/scratch/users/tbellg/pang/grenenet_reads/Table_S5_sample_collection_sequencing_library.csv"
  "/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/samples_data_fix57.csv"
  "/global/scratch/projects/fc_moilab/projects/grenenet-phase1/drive_zenodo/data-intermediate/bioclimvars_experimental_sites_era5.csv"
  "/global/scratch/users/tbellg/gea_grene-net/ARCHIVE/linages_wza_picmin/kendall_0_w_id_n_blocks.csv"
)
for f in "${EXT[@]}"; do
  echo "    - $f"
  $RCLONE copyto "$f" "$DEST/external/$(basename "$f")" "${FLAGS[@]}"
done
# TAIR10 annotation dir (gene/TE GFFs)
$RCLONE copy "/global/home/users/tbellg/ara_key_files" "$DEST/external/ara_key_files" "${FLAGS[@]}"

# 3. The GEA results (the big one: af_store 52G + matrices + outputs = ~67G)
echo "### [3/3] results/grenenet_gea -> results_grenenet_gea/  (~67 GB)"
$RCLONE copy "$PROJ/results/grenenet_gea" "$DEST/results_grenenet_gea" "${FLAGS[@]}"

echo "### finished: $(date)"
echo "### verify:  rclone size \"$DEST\""
