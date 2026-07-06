#!/usr/bin/env python3
"""Discriminate the two explanations for "SV-enriched in JOINT (any-garden selection) but
flat across CLIMATE": genuine per-garden local adaptation vs a directionless per-garden
z-variance inflation in SV-containing blocks (both predict JOINT-up / CLIMATE-flat).

Reuses Z/C from multisite_founder_gwas_clq90_pc1.npz (MAF>=1%) + sv_landscape_clq0.9.csv.
Three reads:
  (1) Per-garden SV z-variance excess (size-matched): uniform ~1.03 across all 30 gardens =>
      systematic class effect (variance-inflation-like); spiky at a few gardens => site-
      specific selection (biology-like). Also its correlation with the climate axis (gradient)
      and |axis| (extremes) -- does the per-garden SV excess have ANY climate structure?
  (2) Concentration: for top-JOINT SV blocks, what share of chi2_joint sits in the single
      largest garden? A locus really selected at one place concentrates; diffuse variance
      inflation spreads ~evenly (1/30). Compared to size-matched SNP blocks + a noise sim.
  (3) Locality-index null: is the genome-wide median locality (~0.95) actually elevated over
      the pure-noise expectation (residual df / joint df = 28/30), or right at it?

Run in `basic` env from the repo root. Deterministic (seed=0).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
import lib

OUT = Path("results/grenenet_gea/sv_adaptive")
NPERM = 5000
EDGES = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 70, 100, 150, 250, 10**9]
rng = np.random.default_rng(0)

raw = lib.maf_filter(lib.multisite_gwas_raw("clq90_pc1"), 0.01)
Z, C, Cinv, sites, bio1 = raw["Z"], raw["C"], raw["Cinv"], raw["sites"], raw["bio1"]
S = raw["S"]
# whiten z per garden so "variance" is on the same footing the omnibus uses:
# JOINT = z C^-1 z; the per-garden marginal that JOINT sums is the whitened coordinate.
# Use the raw per-site z^2 (already per-site lambda~1 by construction) for an interpretable
# per-garden read, and chi2_joint for the omnibus.
chi2_joint = raw["chi2_joint"]

# ---- aggregate markers -> block: max|z| per garden, max chi2_joint (matches persite script) ----
absZ = np.abs(Z)
dfz = pd.DataFrame(absZ, columns=[f"s{int(s)}" for s in sites])
dfz["unit"] = raw["unit"]; dfz["chi2_joint"] = chi2_joint
agg = dfz.groupby("unit").max()

L = pd.read_csv(OUT / "sv_landscape_clq0.9.csv")
B = agg.merge(L[["block_id", "n_kept", "has_sv"]], left_index=True, right_on="block_id")
B["bin"] = pd.cut(B.n_kept, EDGES, right=False, labels=False)
B = B[B.bin.notna()].reset_index(drop=True); B["bin"] = B["bin"].astype(int)
bins = B["bin"].to_numpy()
sv_idx = np.where(B.has_sv.to_numpy() == 1)[0]
print(f"{len(B):,} blocks | {len(sv_idx):,} SV-bearing | garden count S={S}")

# ---- (1) per-garden size-matched SV z^2 excess ----
print("\n(1) per-garden SV z^2 excess (size-matched):  ratio = SV / matched-SNP  (>1 = SV noisier/selected there)")
rows = []
for s in sites:
    z2 = (B[f"s{int(s)}"].to_numpy()) ** 2                # (max|z| per block at garden s)^2
    obs, null, ratio, p = lib.matched_perm_test(z2, bins, sv_idx, NPERM)
    rows.append(dict(site=int(s), axis=float(bio1[list(sites).index(s)]), ratio=ratio, p=p))
R = pd.DataFrame(rows)
R["abs_axis"] = R.axis.abs()
print(f"  ratios: median {R.ratio.median():.3f}, range [{R.ratio.min():.3f}, {R.ratio.max():.3f}], "
      f">1 in {int((R.ratio>1).sum())}/{S}, sig(p<.05) in {int((R.p<0.05).sum())}/{S}")
print(f"  spread across gardens (SD of ratio) = {R.ratio.std():.3f}  "
      f"(small => uniform/systematic; large => spiky/site-specific)")
rg = stats.spearmanr(R.axis, R.ratio); ra = stats.spearmanr(R.abs_axis, R.ratio)
print(f"  ratio vs climate axis (gradient): Spearman rho={rg.correlation:+.3f} p={rg.pvalue:.3f}")
print(f"  ratio vs |climate axis| (extremes): Spearman rho={ra.correlation:+.3f} p={ra.pvalue:.3f}")
R.sort_values("axis").to_csv(OUT / "sv_joint_persite_variance.csv", index=False)

# ---- (2) concentration: share of a block's chi2_joint in its single largest garden ----
# per-marker whitened per-garden contribution: use W = C^-1/2 z so sum of squares = chi2_joint
evals, evecs = np.linalg.eigh(Cinv)
Chalf = evecs @ np.diag(np.sqrt(np.clip(evals, 0, None))) @ evecs.T   # C^-1/2 (symmetric)
W = Z @ Chalf                                                        # (M, S), rowsum sq = chi2_joint
w2 = W ** 2
top_share_marker = w2.max(1) / np.clip(w2.sum(1), 1e-12, None)       # per marker
# aggregate to block via the lead JOINT marker per unit
lead = pd.DataFrame({"unit": raw["unit"], "chi2": chi2_joint, "share": top_share_marker})
lead = lead.loc[lead.groupby("unit")["chi2"].idxmax()].set_index("unit")
B2 = B.merge(lead[["share"]], left_on="block_id", right_index=True)
# noise expectation: max of S iid chi2_1 over their sum
sim = rng.chisquare(1, size=(20000, S)); noise_share = (sim.max(1) / sim.sum(1)).mean()
top_by_joint = B2.nlargest(int(0.02 * len(B2)), "chi2_joint")
sv_share = top_by_joint.loc[top_by_joint.has_sv == 1, "share"].mean()
snp_share = top_by_joint.loc[top_by_joint.has_sv == 0, "share"].mean()
print(f"\n(2) top-2% JOINT blocks -- share of chi2_joint in the single largest garden:")
print(f"  SV blocks {sv_share:.3f} | SNP blocks {snp_share:.3f} | pure-noise expectation {noise_share:.3f}")
print(f"  (higher than noise => concentrated at one garden = locus-like; ~noise => diffuse)")

# ---- (3) locality-index null: observed median vs pure-noise (z ~ MVN(0,C)) ----
loc_obs = np.nanmedian(np.clip(chi2_joint - raw["z_global"]**2 - raw["z_clim"]**2, 0, None) / chi2_joint)
Lchol = np.linalg.cholesky(C + 1e-9*np.eye(S))
one = raw["one"]; dg = raw["dg"]
c0 = (bio1 - bio1.mean())/bio1.std(); cc = c0 - (float(one@Cinv@c0)/dg)*one
zn = (Lchol @ rng.standard_normal((S, 40000))).T                    # null z ~ MVN(0,C)
gj = np.einsum("mi,ij,mj->m", zn, Cinv, zn)
gg = (zn @ Cinv @ one)/np.sqrt(dg); gcl = (zn @ Cinv @ cc)/np.sqrt(float(cc@Cinv@cc))
loc_null = np.median(np.clip(gj - gg**2 - gcl**2, 0, None)/gj)
print(f"\n(3) locality index: observed median {loc_obs:.4f} | pure-noise median {loc_null:.4f} "
      f"| naive df ratio {(S-2)/S:.4f}")
print(f"  => observed is {'ELEVATED over' if loc_obs > loc_null + 0.005 else 'AT'} the noise floor "
      f"({'excess local signal genome-wide' if loc_obs > loc_null + 0.005 else 'residual ~ pure noise for the bulk'})")

print("\n[done] -> sv_joint_persite_variance.csv")
