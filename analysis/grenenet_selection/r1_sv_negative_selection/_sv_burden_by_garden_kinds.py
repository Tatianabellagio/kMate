#!/usr/bin/env python
"""Per-garden burden LMM for EVERY variant kind, not just the headline CDS-insertion model.

`_sv_burden_lmm.py` fits the headline model (shared CDS insertions) per garden and writes
sv_burden_lmm_by_garden{TAG}.csv. This script runs the identical model -- same design matrix, same
covariates, same kinship random effect, same LMM class, imported not copied -- for insertions,
deletions, short indels and SNPs, so the per-garden effects can be compared across kinds:

  y_f(garden s) = log((h0_f + dh_sf + eps)/(h0_f + eps))
  y ~ 1 + z(rate kind_shared_CDS_all) + z(rate kind_shared_intergenic_all)
        + z(callrate_sv) + z(callrate_snp) + cactus + u,   u ~ N(0, s2g K)

beta = the CDS term, so it is a CDS-vs-intergenic contrast of the same kind and origin class, with
relatedness controlled. A CONTEXT env var overrides the tested context (default CDS).

Env: kmate (needs SV_UNIT=event for the Truvari-collapsed event unit, matching the notebooks).
Writes results/sv_adaptive/sv_burden_lmm_by_garden_kinds{TAG}.csv
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _sv_burden_lmm as BL          # LMM class, OUT, EF, EPS, TAG -- no side effects on import

CONTEXT = os.environ.get("CONTEXT", "CDS")
KINDS = ("ins", "del", "indel", "snp")


def main():
    D = pd.read_csv(f"{BL.OUT}/founder_burden{BL.TAG}.csv", dtype={"founder": str})
    K = np.load(f"{BL.OUT}/founder_kinship.npy")
    Zf = np.load(BL.EF, allow_pickle=True)
    assert np.array_equal(Zf["founders"].astype(str), D.founder.to_numpy())
    Y = np.log((Zf["H0"][None, :] + Zf["DH_rel"] + BL.EPS) / (Zf["H0"][None, :] + BL.EPS))
    L = BL.LMM(K)
    cov = np.column_stack([np.ones(len(D)), BL.z(D.callrate_sv), BL.z(D.callrate_snp), D.is_cactus.astype(float)])
    rate = lambda c: (D[c] / D[f"called_{c}"]).to_numpy(float)

    rows = []
    for kind in KINDS:
        t = rate(f"{kind}_shared_{CONTEXT}_all"); ref = rate(f"{kind}_shared_intergenic_all")
        if t.std() == 0 or ref.std() == 0:
            print(f"[skip] {kind}: no variance"); continue
        X = np.column_stack([cov, BL.z(t), BL.z(ref)])
        b_all, se_all, _ = L.fit(Y.mean(0), X)
        for s in range(Y.shape[0]):
            b, se, h2 = L.fit(Y[s], X)
            rows.append(dict(kind=kind, context=CONTEXT, site=int(Zf["sites"][s]), bio1=Zf["bio1"][s],
                             zone=str(Zf["zone"][s]), beta=b[-2], se=se[-2],
                             p=2 * stats.norm.sf(abs(b[-2] / se[-2])), h2=h2,
                             beta_meanfitness=b_all[-2], se_meanfitness=se_all[-2]))
        print(f"[{kind}] mean-fitness beta {b_all[-2]:+.4f} +- {se_all[-2]:.4f}", flush=True)
    R = pd.DataFrame(rows)
    f = f"{BL.OUT}/sv_burden_lmm_by_garden_kinds{BL.TAG}.csv"
    R.to_csv(f, index=False)

    print(f"\nper-garden {CONTEXT}-vs-intergenic burden effect, inverse-variance pooled over gardens:")
    for kind, g in R.groupby("kind", sort=False):
        b = g.beta.to_numpy(); v = g.se.to_numpy() ** 2; w = 1 / v
        pooled = (w * b).sum() / w.sum(); pse = np.sqrt(1 / w.sum())
        Q = (w * (b - pooled) ** 2).sum(); df = len(b) - 1
        I2 = max(0.0, (Q - df) / Q) * 100
        rho = stats.spearmanr(g.bio1, g.beta)
        print(f"  {kind:<6} pooled {pooled:+.4f} +- {pse:.4f} (p = {2*stats.norm.sf(abs(pooled/pse)):.2g}); "
              f"negative in {(b<0).sum()}/{len(b)}; Q = {Q:.1f}/{df} df (p = {stats.chi2.sf(Q,df):.2f}), "
              f"I2 = {I2:.0f}%; beta vs bio1 rho = {rho.statistic:+.2f} (p = {rho.pvalue:.2g})")
    print(f"[wrote] {f}")


if __name__ == "__main__":
    main()
