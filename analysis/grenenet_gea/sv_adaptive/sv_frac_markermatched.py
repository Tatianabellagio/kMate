#!/usr/bin/env python3
"""Does the original SV-frac JOINT enrichment survive matching on hap-cluster count?

sv_joint_mechanism.py showed the has_sv<->JOINT association is a per-block MAX-over-hap-
cluster-markers aggregation effect: block JOINT score = min-p over the block's n_markers
tested hap-clusters, so SV-bearing blocks (more haplotype diversity -> more clusters) get a
lower min-p for free. Matching on n_markers collapsed the has_sv excess (x1.05 -> x0.999).

This re-tests the NOTEBOOK HEADLINE metric -- sv_frac (SV share of a block's records) among
top-JOINT blocks -- three ways on the SAME clq0.9 blocks:
  (A) ORIGINAL: size-matched on n_kept only            (reproduces sv_enrichment.csv)
  (B) + hap-cluster count: matched on n_kept x n_markers (the aggregation confound)
  (C) Šidák-ranked: rank blocks by 1-(1-min_p)^n_markers (conservative per-block multiple-
      testing correction at the source), then the ORIGINAL n_kept matching.
If (B)/(C) collapse the enrichment toward 1, the headline is largely the aggregation artifact;
if it persists, the SV-frac signal is robust to it.

n_markers here = # tested hap-cluster markers per block (rows per unit in the GWAS csv).
Run in `basic` env from the repo root. Deterministic (seed=0). NPERM=10000 (matches original).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pathlib import Path
import numpy as np, pandas as pd
import lib

OUT = Path("results/grenenet_gea/sv_adaptive")
NPERM, TOP = 10000, [0.005, 0.01, 0.02]
EDGES = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 70, 100, 150, 250, 10**9]
MMAX = 5                                   # n_markers bins: 1,2,3,4,5+ (coarse -> populated cells)

g = pd.read_csv(f"{lib.GEA}/hapfreq/multisite_founder_gwas_clq90_pc1.csv")
gc = g.groupby("unit").agg(p_joint=("p_joint", "min")).reset_index()
gc["n_markers"] = g.groupby("unit").size().values
gc["p_sidak"] = 1 - (1 - gc.p_joint) ** gc.n_markers          # (C) per-block multiple-test correction

L = pd.read_csv(OUT / "sv_landscape_clq0.9.csv")
df = gc.merge(L[["block_id", "n_kept", "n_sv", "sv_frac", "has_sv"]],
              left_on="unit", right_on="block_id", how="inner").reset_index(drop=True)
df["kbin"] = pd.cut(df.n_kept, EDGES, right=False, labels=False)
df = df[df.kbin.notna()].reset_index(drop=True); df["kbin"] = df.kbin.astype(int)
df["mbin"] = np.minimum(df.n_markers, MMAX)
df["jbin"] = df.kbin * 10 + df.mbin                           # joint n_kept x n_markers cell
print(f"{len(df):,} blocks | n_markers median {df.n_markers.median():.0f} mean {df.n_markers.mean():.2f} "
      f"| SV blocks n_markers median {df.loc[df.has_sv==1,'n_markers'].median():.0f}")

kb = df.kbin.to_numpy(); jb = df.jbin.to_numpy()

def run(metric):
    vals = df[metric].to_numpy(float)
    rows = []
    for frac in TOP:
        n = int(round(frac * len(df)))
        idx_j = df.p_joint.to_numpy().argsort()[:n]           # top by JOINT min-p (A,B)
        idx_s = df.p_sidak.to_numpy().argsort()[:n]           # top by Šidák block-p (C)
        # imbalance the confound creates: mean n_markers in the top set vs genome-wide
        nm_top, nm_all = df.n_markers.iloc[idx_j].mean(), df.n_markers.mean()
        # min control-pool size across the top set's joint cells (degeneracy guard)
        from collections import Counter
        cell_n = Counter(jb.tolist()); min_pool = min(cell_n[jb[i]] for i in idx_j)
        oA, nA, rA, pA = lib.matched_perm_test(vals, kb, idx_j, NPERM)
        oB, nB, rB, pB = lib.matched_perm_test(vals, jb, idx_j, NPERM)
        oC, nC, rC, pC = lib.matched_perm_test(vals, kb, idx_s, NPERM)
        rows.append(dict(top=f"{frac:.1%}", n=n, nm_top=round(nm_top, 2), nm_all=round(nm_all, 2),
                         min_pool=min_pool,
                         A_kept=f"x{rA:.2f} p={pA:.3f}", B_kept_x_nmark=f"x{rB:.2f} p={pB:.3f}",
                         C_sidak=f"x{rC:.2f} p={pC:.3f}"))
    return pd.DataFrame(rows)

print("\n================  metric = sv_frac (the notebook headline)  ================")
R = run("sv_frac"); print(R.to_string(index=False))
R.to_csv(OUT / "sv_frac_markermatched.csv", index=False)

# validate (A) reproduces the saved original
orig = pd.read_csv(OUT / "sv_enrichment.csv")
orig = orig[(orig.contrast == "JOINT") & (orig.metric == "sv_frac")].set_index("top")
mine = {r.top: float(r.A_kept.split()[0][1:]) for _, r in R.iterrows()}
print("\nvalidation vs saved sv_enrichment.csv (JOINT sv_frac enrich):")
for t in ["0.5%", "1.0%", "2.0%"]:
    print(f"  {t}: mine(A) x{mine[t]:.2f}  saved x{orig.loc[t,'enrich']:.2f}  "
          f"{'OK' if abs(mine[t]-orig.loc[t,'enrich'])<0.06 else 'MISMATCH'}")

print("\n================  metric = has_sv (binary, for reference)  ================")
print(run("has_sv").to_string(index=False))
print(f"\n[done] -> {OUT}/sv_frac_markermatched.csv")
