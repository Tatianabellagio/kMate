#!/usr/bin/env python
"""Audit rdavp_varpart.py: verify the RDA constrained-variance + permutation math
with positive/negative controls and identity checks, so we trust the null result."""
import numpy as np, pandas as pd, sys
sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel")
import rdavp_varpart as R
n = R.n
rng = np.random.RandomState(1)

print(f"n sites = {n}")

# ---- Test 1: Gram-trick identity  trace(Hx G) == ||Hx Yz||^2_F ----
Y = rng.randn(n, 500); X = rng.randn(n, 3); Z = rng.randn(n, 3)
G, tot, Hz = R._prep(Y, Z)
Xz = X - Hz @ X
Hx = Xz @ np.linalg.pinv(Xz.T @ Xz) @ Xz.T
Yz = (Y - Y.mean(0)) - Hz @ (Y - Y.mean(0))
direct = np.sum((Hx @ Yz) ** 2)              # ||Hx Yz||^2_F
trick = R._con(X, Hz, G)                      # trace(Hx G)
print(f"\n[1] Gram-trick identity: direct={direct:.4f} trick={trick:.4f} "
      f"-> {'OK' if abs(direct-trick)<1e-6 else 'FAIL'}")
print(f"    total SS: trace(G)={tot:.4f}  ||Yz||^2={np.sum(Yz**2):.4f} "
      f"-> {'OK' if abs(tot-np.sum(Yz**2))<1e-6 else 'FAIL'}")

# ---- Test 2: POSITIVE control -- Y genuinely driven by climate PC1 ----
sp = R.pools.drop_duplicates("site").set_index("site").loc[R.sites]
B = sp[[b for b in R.BIOS if b in sp.columns]].to_numpy(float); Bz=(B-B.mean(0))/B.std(0)
Ub, Sb, _ = np.linalg.svd(Bz - Bz.mean(0), full_matrices=False)
CLIM = Ub[:, :3] * Sb[:3]
sig = CLIM[:, 0:1] @ rng.randn(1, 400) + 0.3 * rng.randn(n, 400)   # 400 loci ~ climate PC1
r2_pos, p_pos = R.rda(sig, CLIM, None, perm=True)
print(f"\n[2] POSITIVE control (Y = climate-PC1 + noise): R2={r2_pos:.3f} perm_p={p_pos:.4f} "
      f"-> {'OK (detects signal)' if p_pos<0.01 and r2_pos>0.5 else 'FAIL'}")

# ---- Test 3: NEGATIVE control -- Y random, unrelated to climate ----
r2s, ps = [], []
for seed in range(5):
    rr = np.random.RandomState(seed)
    Yn = rr.randn(n, 2000)
    r2, p = R.rda(Yn, CLIM, None, perm=True)
    r2s.append(r2); ps.append(p)
exp = 3.0 / (n - 1)   # E[R2] under null with p=3 predictors
print(f"\n[3] NEGATIVE control (Y random, 5 seeds): mean R2={np.mean(r2s):.3f} "
      f"(null expectation p/(n-1)={exp:.3f})  perm_p range=[{min(ps):.2f},{max(ps):.2f}] "
      f"-> {'OK (R2~null, p~uniform)' if abs(np.mean(r2s)-exp)<0.03 and min(ps)>0.1 else 'CHECK'}")

# ---- Test 4: site-AF aggregation cross-check (SV, small) ----
af = np.load(f"{R.CM}/sv_gen9_af.npy")
S = R.site_af("sv")
# manual for first site, first 3 variants
s0 = R.sites[0]; m = R.site_of == s0; ws = R.w[m] / R.w[m].sum()
manual = ws @ af[m, :3].astype(np.float64)
print(f"\n[4] site-AF aggregation (site {s0}, 3 variants): script={S[0,:3]}  manual={manual} "
      f"-> {'OK' if np.allclose(S[0,:3], manual) else 'FAIL'}")

# ---- Test 5: partition sanity -- full R2 >= pure R2 (confounded >= 0) ----
Zst = None
Ssnp = S  # reuse sv as a quick response
r2_full, _ = R.rda(Ssnp, CLIM, None, perm=False)
r2_pure, _ = R.rda(Ssnp, CLIM, CLIM[:, :2] + rng.randn(n, 2) * 0.1, perm=False)  # arbitrary Z
print(f"\n[5] partition sanity: full R2={r2_full:.3f} >= pure(cond) R2={r2_pure:.3f} ? "
      f"-> {'OK' if r2_full >= r2_pure - 1e-9 else 'note: depends on Z'}")

print("\nAUDIT DONE")
