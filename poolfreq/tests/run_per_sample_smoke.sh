#!/bin/bash
#SBATCH --job-name=per_sample_smoke
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=logs/per_sample_smoke_%j.out
#SBATCH --error=logs/per_sample_smoke_%j.err

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate/poolfreq
mkdir -p results/smoke

PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

# Run per-sample driver on skewed5 simulation (clean known truth)
echo "[$(date)] === skewed5 simulation ==="
$PYTHON src/per_sample_driver.py \
    --cn-kmer-prefix data/cn_full \
    --cn-var data/cn_var_82.cn_var.npz \
    --cn-var-meta data/cn_var_82.meta.npz \
    --reads data/sim_chr1_skewed/skewed5_pool_1.fq.gz data/sim_chr1_skewed/skewed5_pool_2.fq.gz \
    --sample skewed5_sim \
    --out results/smoke/skewed5_sim.tsv \
    --threads 8

# Run on SEEDMIX_S1 (real data)
echo ""
echo "[$(date)] === SEEDMIX_S1 ==="
$PYTHON src/per_sample_driver.py \
    --cn-kmer-prefix data/cn_full \
    --cn-var data/cn_var_82.cn_var.npz \
    --cn-var-meta data/cn_var_82.meta.npz \
    --reads /global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz \
            /global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz \
    --sample SEEDMIX_S1 \
    --out results/smoke/SEEDMIX_S1.tsv \
    --threads 8

echo ""
echo "[$(date)] DONE"
ls -la results/smoke/

# Quick analysis
$PYTHON <<'PYEOF'
import pandas as pd
import numpy as np
print("\n=== Output validation ===")
for sample in ["skewed5_sim", "SEEDMIX_S1"]:
    df = pd.read_csv(f"results/smoke/{sample}.tsv", sep="\t")
    print(f"\n{sample}:")
    print(f"  records: {len(df):,}")
    print(f"  alt_freq: min={df.alt_freq.min():.4f}  max={df.alt_freq.max():.4f}  "
          f"mean={df.alt_freq.mean():.4f}  median={df.alt_freq.median():.4f}")
    print(f"  records with alt_freq > 0.01: {(df.alt_freq > 0.01).sum():,}")
    print(f"  records with alt_freq > 0.1:  {(df.alt_freq > 0.1).sum():,}")
PYEOF
