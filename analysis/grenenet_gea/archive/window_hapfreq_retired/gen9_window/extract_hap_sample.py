#!/usr/bin/env python
"""Collapse one sample's window founder-h into the 38,871 haplotype-cluster frequencies.
For trajectory t (unit u, founder set member): freq = sum_{f in member} h[u, f].
Output: hap/persample/{sample}.npy  float32 [n_traj] (traj_id order, chr1..chr5)."""
import os, sys
import numpy as np

ROOT = "/global/scratch/users/tbellg/kmate"
WIN = f"{ROOT}/results/grenenet_kmate_window"
GW = f"{ROOT}/results/grenenet_gea/gen9_window"
sample = sys.argv[1]
out = f"{GW}/hap/persample/{sample}.npy"
if os.path.exists(out):
    print(f"{sample}: exists"); sys.exit(0)

vecs = []
for ci, chrlc in enumerate(["chr1", "chr2", "chr3", "chr4", "chr5"], 1):
    npz = np.load(f"{WIN}/{sample}_Chr{ci}.h_blocks_per_chrom.npz")
    h = npz[f"Chr{ci}_h_blocks"]                       # [n_units x 231]
    mem = np.load(f"{GW}/hap/{chrlc}_member.npz")
    member = mem["member"]; uidx = mem["unit_idx"]     # [n_traj x 231], [n_traj]
    Hexp = h[uidx]                                      # [n_traj x 231]
    vecs.append((Hexp * member).sum(1).astype(np.float32))
v = np.concatenate(vecs)
os.makedirs(f"{GW}/hap/persample", exist_ok=True)
tmp = out + ".tmp.npy"; np.save(tmp, v); os.replace(tmp, out)
print(f"{sample}: {v.size:,} hap freqs (mean {np.nanmean(v):.4f})")
