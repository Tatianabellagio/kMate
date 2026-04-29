#!/bin/bash
#SBATCH --job-name=per_sample_smoke
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/per_sample_smoke_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/tests/logs/per_sample_smoke_%j.err

set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq
mkdir -p results/smoke

PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python

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
    --reads /home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz \
            /home/tbellagio/scratch/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz \
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
