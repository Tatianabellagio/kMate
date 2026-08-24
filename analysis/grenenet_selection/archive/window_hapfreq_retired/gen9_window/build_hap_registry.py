#!/usr/bin/env python
"""Define the haplotype-trajectory test units from the per-unit cluster registry, and
the founder-membership needed to collapse each sample's window h into cluster frequencies.

A trajectory = a kept haplotype cluster: founding cluster_freq >= 2/231 (>=2 founders) AND
n_sig >= 2 (>=2 signature variants define it). is_major flags the unit's largest cluster
(the natural reference category). hap-cluster frequency for a sample = sum of window h over
the cluster's founders.

Outputs (results/grenenet_gea/gen9_window/hap/):
  registry.csv       traj_id,chrom,unit_idx,start,end,n_eff,cluster,founding_freq,n_sig,is_major
  {chrlc}_member.npz  member [n_traj_chrom x 231] bool, unit_idx [n_traj_chrom], traj_id [n_traj_chrom]
"""
import os, glob
import numpy as np, pandas as pd

GW = "results/grenenet_gea/gen9_window"
OUT = f"{GW}/hap"
os.makedirs(OUT, exist_ok=True)
MINFREQ, MINSIG = 2 / 231, 2

reg = []
tid = 0
for chrlc in ["chr1", "chr2", "chr3", "chr4", "chr5"]:
    C = pd.read_csv(f"{GW}/clusters/{chrlc}_clusters.csv")
    lab = np.load(f"{GW}/clusters/{chrlc}_labels.npz")
    labels = lab["labels"]                                   # [n_units x 231]
    kept = C[(C.cluster_freq >= MINFREQ) & (C.n_sig >= MINSIG)].copy()
    # is_major = largest founding-freq cluster within the unit (over ALL clusters, kept or not)
    major = C.sort_values("cluster_freq").groupby("unit_idx").tail(1)[["unit_idx", "cluster"]]
    major_set = set(map(tuple, major.values))
    mem, uidx, tids = [], [], []
    rows = []
    for _, r in kept.iterrows():
        ui, cl = int(r.unit_idx), int(r.cluster)
        m = labels[ui] == cl                                 # 231-bool founder membership
        mem.append(m); uidx.append(ui); tids.append(tid)
        rows.append((tid, r.chrom, ui, int(r.start), int(r.end), float(r.n_eff), cl,
                     float(r.cluster_freq), int(r.n_sig), (ui, cl) in major_set))
        tid += 1
    reg.extend(rows)
    np.savez(f"{OUT}/{chrlc}_member.npz",
             member=np.array(mem, bool), unit_idx=np.array(uidx, np.int32),
             traj_id=np.array(tids, np.int32))
    print(f"{chrlc}: {len(kept):,} trajectories", flush=True)

R = pd.DataFrame(reg, columns=["traj_id", "chrom", "unit_idx", "start", "end", "n_eff",
                               "cluster", "founding_freq", "n_sig", "is_major"])
R.to_csv(f"{OUT}/registry.csv", index=False)
print(f"total {len(R):,} hap trajectories ({(~R.is_major).sum():,} non-major) -> {OUT}/registry.csv")
