#!/usr/bin/env python
"""Founding (gen0) haplotype-cluster frequencies from the 8 SEEDMIX window-mode reps.
p0 = mean over reps of (sum window-h over each cluster's founders); v0 = among-rep
variance of logit(p0) / n_reps (the gen0 point variance for Pipeline B)."""
import os, glob, sys
import numpy as np
GW = "results/grenenet_gea/gen9_window"
SM = "results/grenenet_kmate_window_seedmix"
EPS = 1e-3
def logit(p):
    p = np.clip(p, EPS, 1 - EPS); return np.log(p / (1 - p))

reps = sorted(glob.glob(f"{SM}/SEEDMIX_S*.tsv"))
reps = [os.path.basename(f)[:-4] for f in reps if "_Chr" not in os.path.basename(f)]
mats = []
for s in reps:
    vecs = []
    for ci, chrlc in enumerate(["chr1", "chr2", "chr3", "chr4", "chr5"], 1):
        h = np.load(f"{SM}/{s}_Chr{ci}.h_blocks_per_chrom.npz")[f"Chr{ci}_h_blocks"]
        mem = np.load(f"{GW}/hap/{chrlc}_member.npz")
        vecs.append((h[mem["unit_idx"]] * mem["member"]).sum(1).astype(np.float32))
    mats.append(np.concatenate(vecs))
M = np.vstack(mats)                                  # [n_reps x n_traj]
p0 = M.mean(0)
lv = logit(M).var(0, ddof=1) / M.shape[0]
np.save(f"{GW}/hap/p0_seedmix.npy", p0.astype(np.float32))
np.save(f"{GW}/hap/v0_seedmix.npy", lv.astype(np.float32))
print(f"{M.shape[0]} SEEDMIX reps -> p0 (mean {p0.mean():.4f}) + v0 (median {np.median(lv):.3f}); {p0.size:,} trajectories")
