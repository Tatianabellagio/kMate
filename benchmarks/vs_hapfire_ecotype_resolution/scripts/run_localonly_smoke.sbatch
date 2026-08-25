#!/bin/bash
#SBATCH --job-name=localonly_smoke
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH --output=logs/localonly_smoke_%j.out
#SBATCH --error=logs/localonly_smoke_%j.out
# Smoke-test the new --local-only (global-free) window mode on a g3 recombinant
# p80 pool. Confirms: runs end-to-end, low/empty windows + unassigned records
# emit NaN AF (NOT global_h), HMM smoothing still applies to local fits.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs analysis/panel_qc/localonly_smoke
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
D=benchmarks/p80/sims/cov10_n50_g3_s42_hotspots_p80_chr1
OUT=analysis/panel_qc/localonly_smoke

run() {  # mode-label  extra-flags...
  local label=$1; shift
  echo "=== [$label] $(date) ==="
  $PY src/per_sample_per_chrom.py \
    --kmer-pa-prefix benchmarks/p80/data/kmer_pa_p80_filt2/kmer_pa \
    --var-pa     benchmarks/p80/data/var_pa_p80.var_pa.npz \
    --var-called benchmarks/p80/data/var_pa_p80.var_called.npz \
    --var-meta   benchmarks/p80/data/var_pa_p80.meta.npz \
    --reads $D/reads/r1.fq $D/reads/r2.fq \
    --sample g3_$label --out $OUT/g3_$label.tsv \
    --threads 8 --chroms Chr1 --kmer-weight inv_mb \
    --block-mode window --window-bp 10000 --min-kmers-per-block 200 "$@"
}

# Global-free (new) vs production-window (anchored fallback) for comparison.
run localonly --local-only
run anchored

echo "=== NaN summary (local-only should have NaN AF in low/empty windows) ==="
$PY - <<'PYEOF'
import pandas as pd, numpy as np
for lab in ("localonly","anchored"):
    df = pd.read_csv(f"analysis/panel_qc/localonly_smoke/g3_{lab}.tsv", sep="\t")
    n = len(df); nan = df["alt_freq"].isna().sum()
    print(f"{lab:10s}  records={n:,}  NaN_AF={nan:,} ({100*nan/n:.2f}%)  "
          f"finite_AF_mean={np.nanmean(df['alt_freq']):.4f}")
PYEOF
echo "=== done $(date) ==="
