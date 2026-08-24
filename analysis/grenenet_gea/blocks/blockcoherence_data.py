#!/usr/bin/env python
"""Compute per-block PC1 variance-explained for each CLQcut block map (Chr1),
save a tidy CSV for the coherence notebook. Run in the kmate env.
"""
import numpy as np, pandas as pd
from eval_block_coherence import load_allclass_af, pc1_var_explained

CHR = "Chr1"
BR = "analysis/grenenet_gea/blocks_recompute"
MAPS = {0.5: f"{BR}/chr1_clq0.5_blocks_clq0.5.tsv",
        0.7: f"{BR}/chr1_clq0.7_blocks_clq0.7.tsv",
        0.9: f"{BR}/chr1_clq0.9_blocks_clq0.9.tsv"}

print("loading all-class gen9 pool AF ...")
pos, AF = load_allclass_af(CHR, "gen9")
print(f"  {AF.shape[0]} pools x {AF.shape[1]} records")

rows = []
for clq, bf in MAPS.items():
    bl = pd.read_csv(bf, sep="\t")
    for s, e in zip(bl["start_pos"], bl["end_pos"]):
        lo, hi = np.searchsorted(pos, s), np.searchsorted(pos, e, side="right")
        n = hi - lo
        ve = pc1_var_explained(AF[:, lo:hi]) if n >= 2 else np.nan
        rows.append((clq, s, e, n, ve))
    print(f"  CLQcut={clq}: {len(bl)} blocks done")

df = pd.DataFrame(rows, columns=["clqcut", "start_pos", "end_pos", "n_variants", "pc1_ve"])
out = f"{BR}/chr1_pc1ve_by_clqcut.csv"
df.to_csv(out, index=False)
print(f"wrote {len(df)} block rows -> {out}")
print(df.groupby("clqcut").agg(nblk=("pc1_ve","size"),
      medVE=("pc1_ve","median"), frac_ge0p7=("pc1_ve", lambda x:(x>=0.7).mean())).round(3))
