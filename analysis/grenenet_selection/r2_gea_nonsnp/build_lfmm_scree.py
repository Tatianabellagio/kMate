#!/usr/bin/env python
"""One-time: PCA scree + broken-stick on the SNP Δp matrix -> small CSV.

Loads the 5.6 GB Δp memmap ONCE here so the K-selection notebook can stay
lightweight (it just reads the small CSV this writes). Re-run only if the SNP
Δp matrix changes.
"""
import numpy as np, pandas as pd

ROOT = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection"
STEM = f"{ROOT}/phase1_replication/lfmm/lfmm_snp_gen9"
OUT  = f"{ROOT}/lfmm/scree_snp_gen9.csv"

np_pools, nr = np.loadtxt(f"{STEM}_dims.txt", dtype=int)
np_pools, nr = int(np_pools), int(nr)
# SEQUENTIAL read of the whole 5.6 GB file (fast on shared scratch); reshape in RAM.
# (Scattered column indexing on the memmap is ~71M random reads and crawls on scratch.)
print(f"reading {np_pools}x{nr:,} float64 sequentially ...", flush=True)
Y = np.fromfile(f"{STEM}_Y.f64", dtype=np.float64).reshape(np_pools, nr)
print("loaded; computing covariance over all records", flush=True)
Yc = Y - Y.mean(axis=0, keepdims=True)                   # center each record across pools
C  = (Yc @ Yc.T) / Yc.shape[1]                           # [pools x pools] covariance (all records)
evals = np.clip(np.linalg.eigvalsh(C)[::-1], 0, None)    # descending, non-negative
prop  = evals / evals.sum()
cum   = np.cumsum(prop)
p = len(evals)
bstick = np.array([np.sum(1.0 / np.arange(k, p + 1)) for k in range(1, p + 1)]) / p

pd.DataFrame({"component": np.arange(1, p + 1),
              "eig_prop": prop, "broken_stick": bstick, "cum_var": cum}).to_csv(OUT, index=False)

above = prop > bstick
K_bstick = int(np.argmax(~above)) if (~above).any() else p
print(f"pools={np_pools} records={nr:,} (all records used)")
print(f"broken-stick K (first crossover) = {K_bstick}")
print(f"components above null (total)    = {int(above.sum())}")
print(f"cum var at K=16                  = {cum[15]:.3f}")
print(f"eig_prop K12 -> K16 -> K20: {prop[11]:.4f} -> {prop[15]:.4f} -> {prop[19]:.4f}")
print("wrote", OUT)
