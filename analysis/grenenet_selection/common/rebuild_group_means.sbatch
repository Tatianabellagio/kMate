#!/bin/bash
#SBATCH --job-name=rebuild_gm
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/rebuild_group_means_%j.out
#SBATCH --error=logs/rebuild_group_means_%j.out
# Regenerate the GEA group_means cache from the CURRENT segregating-panel TSVs.
# The old cache (group_means.npz/.pkl, 2026-06-02) was built on the pre-segregating
# 10.33M panel and carried AC=0 monomorphic scaffolding; the per-sample TSVs were
# rebuilt segregating (8.49M) on 06-03. build_group_means now also asserts panel
# alignment + drops any all-absent record (should be 0 on the clean panel).
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
source ~/miniforge3/etc/profile.d/conda.sh
conda activate kmate

echo "host: $(hostname)  start: $(date)"
# remove the stale caches so build_group_means regenerates (it only writes .npz;
# the orphaned .pkl is from an older code path and is removed too)
rm -fv analysis/grenenet_selection/group_means.npz analysis/grenenet_selection/group_means.pkl

cd analysis/grenenet_selection
python - <<'PY'
import lib, numpy as np
df = lib.build_group_means()          # reads segregating TSVs, asserts alignment, caches .npz
print(f"rebuilt group_means: {len(df):,} records, {df.shape[1]} cols")
# sanity: CAM5 window should now be segregating-only (~216, no 302 monomorphic)
w = df[(df.chrom=="Chr2") & (df.pos.between(11531504,11534833))]
print(f"CAM5 gene+/-500 window: {len(w):,} records (expect ~216 segregating, NOT 302)")
PY
echo "done: $(date)"
