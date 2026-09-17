#!/usr/bin/env python
"""Are the EPFL5 and ELIP1 per-garden GWAS hits (Chr3 8.074 / 8.084 Mb) one signal?

Founder LD (r2 and carrier overlap over the 231-founder panel, same order and presence
matrix as the GEMMA scan) between:
  EPFL5_gwas  Chr3:8074451 1->2   (significant in gardens 10,12,27,42)
  EPFL5_gea   Chr3:8074451 1->11  (the GEA record at the same position)
  ELIP1_gwas  Chr3:8084255 1->2   (significant in garden 27)
plus the strongest-LD partners of each in +-30 kb, and the per-garden correlation of their
gen1-3 pool AF change. env kmate, compute node.
"""
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                      # noqa: E402
sys.path.insert(0, f"{lib.GEA}/r3_persite_gwas")
from class_split_gwas import _panel_order, _load_chrom_212      # noqa: E402

VARS = {"EPFL5_gwas": (8074451, 1, 2), "EPFL5_gea": (8074451, 1, 11),
        "ELIP1_gwas": (8084255, 1, 2)}

order, _, _, _ = _panel_order()
pos, rl, al, vp, vc, n_alt, n_cal = _load_chrom_212("chr3", order)
col = {}
for k, (p, r, a) in VARS.items():
    i = np.where((pos == p) & (rl == r) & (al == a))[0]
    print(k, "panel records:", len(i))
    # Chr3:8074451 1->2 has TWO panel records (different ALT bases). Only the second
    # (5 carriers) is the GWAS hit; the first (16 carriers) is significant nowhere.
    col[k] = int(i[-1] if k == "EPFL5_gwas" else i[0])
G = {k: vp[:, c].toarray().ravel().astype(float) for k, c in col.items()}
C = {k: vc[:, c].toarray().ravel().astype(bool) for k, c in col.items()}

print("\nfounders:", len(order))
for k in G:
    print(f"  {k:11s} carriers {int(G[k][C[k]].sum()):3d} / called {int(C[k].sum())}")
print("\npairwise (founders called at both):")
ks = list(G)
for i in range(len(ks)):
    for j in range(i + 1, len(ks)):
        a, b = ks[i], ks[j]
        m = C[a] & C[b]
        x, y = G[a][m], G[b][m]
        r = np.corrcoef(x, y)[0, 1] if x.std() and y.std() else np.nan
        both = int(((x == 1) & (y == 1)).sum())
        print(f"  {a} x {b}: r2={r**2:.3f} (r={r:+.2f})  carriers both={both}, "
              f"only {a}={int(((x==1)&(y==0)).sum())}, only {b}={int(((x==0)&(y==1)).sum())}")

# top LD partners in the window, testable markers only
win = np.where((pos > 8044000) & (pos < 8115000) & (np.minimum(n_alt, len(order) - n_alt) >= 5))[0]
D = vp[:, win].toarray().astype(float)
for k in ["EPFL5_gwas", "ELIP1_gwas"]:
    g = G[k]
    Dz = D - D.mean(0); gz = g - g.mean()
    r = (Dz * gz[:, None]).sum(0) / np.sqrt((Dz ** 2).sum(0) * (gz ** 2).sum() + 1e-12)
    top = np.argsort(-r ** 2)[:8]
    print(f"\ntop LD partners of {k} in 8.044-8.115 Mb:")
    for t in top:
        w = win[t]
        print(f"   {pos[w]:>9d} {rl[w]}->{al[w]}  r2={r[t]**2:.2f}")

# per-garden AF change agreement (gen9 garden means from the pool store)
from build_ownaxis_shortlist import garden_means                # noqa: E402
from screen_3criteria import founding_p0                        # noqa: E402
idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
rows = {}
for k, (p, r, a) in VARS.items():
    s = np.where((idx["chrom"] == "Chr3") & (idx["pos"] == p) & (idx["ref_len"] == r)
                 & (idx["alt_len"] == a))[0]
    rows[k] = int(s[0])
gm = garden_means(np.array(list(rows.values())))
p0 = founding_p0(np.array(list(rows.values())))
dlt = pd.DataFrame({k: gm[s] - p0[i] for i, (k, s) in enumerate(rows.items())})
print("\nper-garden delta (gen9 mean - p0), Spearman:")
print(dlt.corr(method="spearman").round(2))
print(dlt.loc[[10, 12, 27, 42, 55, 4]].round(3))
