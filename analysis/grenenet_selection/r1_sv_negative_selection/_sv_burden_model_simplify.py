#!/usr/bin/env python
"""Which terms of the founder burden model are needed? Coefficients of every covariate, and the CDS
burden effect under a ladder of simpler models.

Outcome: founder mean log-fitness over the 31 gardens (as `_sv_burden_lmm.py`). Full model:
  y ~ 1 + z(CDS rate) + z(intergenic rate) + z(callrate_sv) + z(callrate_snp) + cactus + u,  u ~ N(0, s2g K)
All continuous covariates standardised (per SD); cactus is 0/1 (effect = cactus minus short-read founders).
Every model is fitted twice: with the kinship random effect (REML, the `LMM` class of `_sv_burden_lmm.py`)
and by OLS (no kinship), so the contribution of K is visible too. Models are compared by ML log-likelihood
/ AIC (the LMM refitted by ML for that, since REML likelihoods are not comparable across fixed effects).

Env: kmate, SV_UNIT=event. KINDS (default ins,del,indel,snp).
Writes results/sv_adaptive/sv_burden_model_simplify{TAG}.csv (all coefficients, all models, all kinds).
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats, optimize

os.environ.setdefault("SV_UNIT", "event")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _sv_burden_lmm as BL

KINDS = os.environ.get("KINDS", "ins,del,indel,snp").split(",")
OUTF = f"{BL.OUT}/sv_burden_model_simplify{BL.TAG}.csv"

# model ladder: name -> covariate terms beyond the intercept and the CDS rate
LADDER = [("full", ["intergenic", "callrate_sv", "callrate_snp", "cactus"]),
          ("- callrate_snp", ["intergenic", "callrate_sv", "cactus"]),
          ("- cactus", ["intergenic", "callrate_sv", "callrate_snp"]),
          ("- both call rates", ["intergenic", "cactus"]),
          ("intergenic + callrate_sv", ["intergenic", "callrate_sv"]),
          ("intergenic only", ["intergenic"]),
          ("CDS only", [])]


def ml_fit(L, y, X):
    """ML (not REML) fit of the same LMM, for log-likelihood / AIC across fixed-effect sets."""
    Uy = L.U.T @ y; UX = L.U.T @ X; n = y.size

    def nll(ld):
        d = np.exp(ld); w = 1 / (L.S + d)
        b = np.linalg.solve(UX.T @ (UX * w[:, None]), UX.T @ (Uy * w)); r = Uy - UX @ b
        s2 = (r * r * w).sum() / n
        return 0.5 * (n * np.log(2 * np.pi * s2) + np.log(L.S + d).sum() + n)
    g = np.linspace(-8, 8, 60); v = [nll(q) for q in g]; g0 = g[int(np.argmin(v))]
    return -optimize.minimize_scalar(nll, bounds=(g0 - 0.3, g0 + 0.3), method="bounded").fun


def ols(y, X):
    b, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ b; n, p = X.shape
    s2 = (r @ r) / (n - p); se = np.sqrt(np.diag(np.linalg.inv(X.T @ X)) * s2)
    ll = -0.5 * n * (np.log(2 * np.pi * (r @ r) / n) + 1)
    return b, se, ll


def main():
    D = pd.read_csv(f"{BL.OUT}/founder_burden{BL.TAG}.csv", dtype={"founder": str})
    K = np.load(f"{BL.OUT}/founder_kinship.npy")
    Zf = np.load(BL.EF, allow_pickle=True)
    assert np.array_equal(Zf["founders"].astype(str), D.founder.to_numpy())
    y = np.log((Zf["H0"][None, :] + Zf["DH_rel"] + BL.EPS) / (Zf["H0"][None, :] + BL.EPS)).mean(0)
    L = BL.LMM(K); z = BL.z
    rate = lambda c: (D[c] / D[f"called_{c}"]).to_numpy(float)
    rows = []
    for kind in KINDS:
        T = {"CDS": z(rate(f"{kind}_shared_CDS_all")), "intergenic": z(rate(f"{kind}_shared_intergenic_all")),
             "callrate_sv": z(D.callrate_sv), "callrate_snp": z(D.callrate_snp), "cactus": D.is_cactus.astype(float).to_numpy()}
        for mname, terms in LADDER:
            names = ["CDS"] + terms; X = np.column_stack([np.ones(len(D))] + [T[t] for t in names])
            for fit in ("LMM (kinship)", "OLS (no kinship)"):
                if fit.startswith("LMM"):
                    b, se, h2 = L.fit(y, X); ll = ml_fit(L, y, X); npar = X.shape[1] + 2
                else:
                    b, se, ll = ols(y, X); h2 = np.nan; npar = X.shape[1] + 1
                for i, t in enumerate(names, 1):
                    rows.append(dict(kind=kind, model=mname, fit=fit, term=t, beta=b[i], se=se[i],
                                     p=2 * stats.norm.sf(abs(b[i] / se[i])), h2=h2, loglik_ml=ll, aic=2 * npar - 2 * ll,
                                     n_terms=len(names)))
        print(f"[{kind}] done", flush=True)
    R = pd.DataFrame(rows); R.to_csv(OUTF, index=False)
    print(f"[wrote] {OUTF}")


if __name__ == "__main__":
    main()
