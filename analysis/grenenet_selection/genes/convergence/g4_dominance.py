#!/usr/bin/env python
"""How much of a candidate's per-garden rise is the garden-4 founder sweep?

For each record in a candidate table: gen9 garden mean - p0 per garden, then
  g4_delta          delta in garden 4
  g4set_share       share of the summed positive delta carried by gardens 4, 43, 45, 32
                    (the gardens the round-1 review found moving with the garden-4 founder)
  n_up_outside      gardens OUTSIDE that set with delta >= max(0.05, p0) (allele at least doubled)
  n_dn_outside      gardens outside the set with delta <= -max(0.05, p0/2)
A high share with n_up_outside <= 1 is the garden-4 picture regardless of the climate rho.

Usage: python g4_dominance.py results/screen_top_loci_ownaxis.csv
Output: <table>_g4dom.csv (target_gene + the four columns). env kmate, compute node.
"""
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from build_ownaxis_shortlist import garden_means          # noqa: E402

G4SET = [4, 43, 45, 32]

T = pd.read_csv(sys.argv[1])
gm = garden_means(np.unique(T.store_row.to_numpy()))
rows = []
for r in T.itertuples():
    d = gm[r.store_row] - r.p0
    pos = d.clip(lower=0)
    inset = d.index.isin(G4SET)
    thr_up, thr_dn = max(0.05, r.p0), max(0.05, r.p0 / 2)
    rows.append({"target_gene": r.target_gene,
                 "g4_delta": round(float(d.get(4, np.nan)), 3),
                 "g4set_share": round(float(pos[inset].sum() / max(pos.sum(), 1e-9)), 2),
                 "n_up_outside": int((d[~inset] >= thr_up).sum()),
                 "n_dn_outside": int((d[~inset] <= -thr_dn).sum())})
out = sys.argv[1].replace(".csv", "_g4dom.csv")
pd.DataFrame(rows).to_csv(out, index=False)
print(pd.DataFrame(rows).describe().round(2)); print("wrote", out)
