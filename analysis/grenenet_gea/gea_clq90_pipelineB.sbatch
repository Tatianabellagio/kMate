#!/bin/bash
#SBATCH --job-name=gea_clq90_pB
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=8:00:00
#SBATCH --output=logs/gea_clq90_pB_%j.out
#SBATCH --error=logs/gea_clq90_pB_%j.out
# Pipeline B (pool-seq TEMPORAL climate GEA) on the clq0.9 (r2>=0.9) LD blocks.
#
# KEY: haploblock frequency is a linear PROJECTION of founder (ecotype) frequency h,
# exactly like the per-variant AF. We project the genome-GLOBAL per-sample founder h
# ({ch}_global_h, already estimated) onto the clq0.9 founder->cluster membership
# (hap_membership/*_hapmemb_clq90.npz, already built). This needs NO window-EM re-run;
# global h is partition-independent so it maps onto any block grid. (window h == global
# here, corr 0.9998 — recombination rare; the standing decision is GLOBAL AF + fine blocks
# as test units.)
#
# Chain (all kmate env):
#   1 project cohort  -> hapfreq_clq90/{chunk, hapfreq_matrix.npy, registry, p0_panel}
#   2 p0/v0 seedmix   -> hapfreq_clq90/hapfreq_p0_seedmix.npy (8 SEEDMIX reps, global h)
#   3 trajectories    -> hapfreq_clq90/pipelineB_varlen/traj_{p,v,sites}.npy  (varlen)
#   4 hap-GEA (Sidak) -> .../block_gea.csv, hap_gea.csv
#   5 block-WZA       -> .../block_wza.csv
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
KPY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
G=analysis/grenenet_gea
HFD=results/grenenet_gea/hapfreq_clq90
export HF_DIR="$HFD" MEMB_TAG=clq90 H_SOURCE=global TRAJ_MODE=varlen
export PB_DIR="$HFD/pipelineB_varlen" N_PERM=5000
mkdir -p "$HFD"

echo "[$(date)] 1/5 project cohort (clq90, global h)"
$KPY -u $G/build_hapfreq_matrix.py --memb-tag clq90 --h-source global \
     --out-dir "$HFD" --start 0 --end 0
$KPY -u $G/build_hapfreq_matrix.py --gather --memb-tag clq90 \
     --chunk-dir "$HFD" --out-dir "$HFD"

echo "[$(date)] 2/5 p0/v0 seedmix (clq90, global h)"
$KPY -u $G/build_hapfreq_p0_seedmix.py

echo "[$(date)] 3/5 haplotype trajectories (varlen)"
$KPY -u $G/build_hap_trajectories.py

echo "[$(date)] 4/5 hap-GEA + block_gea (Sidak, ${N_PERM} perms)"
$KPY -u $G/build_hap_gea.py

echo "[$(date)] 5/5 block-WZA"
$KPY -u $G/build_hap_wza.py

echo "[$(date)] DONE -> $PB_DIR/{block_gea.csv,block_wza.csv,hap_gea.csv}"
