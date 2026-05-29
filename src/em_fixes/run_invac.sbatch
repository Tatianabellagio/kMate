#!/bin/bash
#SBATCH --job-name=invac
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/scratch/h_fixes/invac/invac_%j.log

set -euo pipefail
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
REPO=/global/scratch/users/tbellg/kmate
CNPREF=$REPO/data/cn_full_231_v3qc_v3_filt2/cn_Chr1
OUTDIR=$REPO/scratch/h_fixes/invac
G0=$REPO/scratch/g0_sweep_h_test
SM=$REPO/scratch/seedmix_h_test
mkdir -p $OUTDIR

GAMMAS="0.0,0.25,0.5,0.75,1.0"

# --- g0: 5 sims ---
$PY $REPO/src/em_fixes/invac.py \
  --cn-prefix $CNPREF --gammas $GAMMAS \
  --out-prefix $OUTDIR/g0 \
  --counts \
    $G0/filt2_g0_n231_rep0_rand.counts.npy \
    $G0/filt2_g0_n200_rep0_rand.counts.npy \
    $G0/filt2_g0_n50_rep0_cact.counts.npy \
    $G0/filt2_g0_n50_rep1_bal.counts.npy \
    $G0/filt2_g0_n50_rep2_pg.counts.npy \
  --labels \
    g0_n231_rep0_rand g0_n200_rep0_rand g0_n50_rep0_cact g0_n50_rep1_bal g0_n50_rep2_pg

# --- SEEDMIX: 8 reps ---
$PY $REPO/src/em_fixes/invac.py \
  --cn-prefix $CNPREF --gammas $GAMMAS \
  --out-prefix $OUTDIR/seedmix \
  --counts \
    $SM/filt2_S1.counts.npy $SM/filt2_S2.counts.npy $SM/filt2_S3.counts.npy \
    $SM/filt2_S4.counts.npy $SM/filt2_S5.counts.npy $SM/filt2_S6.counts.npy \
    $SM/filt2_S7.counts.npy $SM/filt2_S8.counts.npy \
  --labels S1 S2 S3 S4 S5 S6 S7 S8

echo "DONE invac sweep"
