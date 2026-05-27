#!/bin/bash
#SBATCH --job-name=subsamp_protect
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/subsamp_protect_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/subsamp_protect_%j.err
set -euo pipefail
BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
IN=$BASE/poolfreq/data/cn_full_231_v3qc_v3_filt2/cn_Chr1
# filt2 bins: [2,3,5,11,26,51,101,201,232] -> bin0=ac2, bin1=ac3-4, bin2=ac5-10,...
# protect1: keep all ac=2 (doubletons), subsample ac>=3 to median
# protect2: keep all ac=2 and ac=3-4, subsample ac>=5 to median
for P in 1 2; do
  OUT=$BASE/poolfreq/data/cn_full_231_v3qc_v3_subsampProtect${P}_refilt2/cn_Chr1
  mkdir -p "$(dirname "$OUT")"
  echo "[$(date)] subsampMedian protect-bins=$P -> $OUT"
  $PY -u $BASE/poolfreq/src/build_subsampled_cn.py \
      --in-cn ${IN}.cn.npz --in-meta ${IN}.meta.npz \
      --out-prefix $OUT --target median --seed 42 \
      --protect-bins $P --refilt2
done
echo "[$(date)] DONE"
