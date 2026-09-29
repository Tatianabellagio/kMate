#!/usr/bin/env python
"""Founder-label permutation null for the per-garden burden betas pooled across gardens.

`sv_climate_maf.ipynb` pools the 31 per-garden burden betas (`sv_burden_lmm_by_garden_kinds{TAG}.csv`)
by inverse variance as if the gardens were independent. They are not: every garden's fit uses the same
231 founders and the same burden vector, and founder log-fitness is correlated between gardens. This
script gives the reference distribution that keeps that dependence.

Null: the CDS-burden residual (after regressing the CDS rate on the intergenic rate of the same kind, as
`_sv_burden_lmm.py` does for its headline permutation) is permuted across founders ONCE per draw and
the same permuted vector is used in all 31 gardens; the 31 fits, the pooled beta, the inverse-variance
meta-regression slope on garden bio1 and the mean-fitness beta are all recomputed on it.
Same model, same LMM class, same covariates as `_sv_burden_by_garden_kinds.py` (imported).

Env: kmate, SV_UNIT=event. NPERM (default 3000), KINDS (default "ins,del").
Writes results/sv_adaptive/sv_burden_pool_null{TAG}.npz (per kind: observed per-garden betas / SEs,
null per-garden betas, null pooled betas, null slopes, null mean-fitness betas) and
sv_burden_pool_null{TAG}.csv (one summary row per kind). Skips a kind already present with >= NPERM
draws unless FORCE=1.
"""
import os, sys, time
import numpy as np, pandas as pd
from scipy import stats

os.environ.setdefault("SV_UNIT", "event")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _sv_burden_lmm as BL

NPERM = int(os.environ.get("NPERM", 3000))
KINDS = os.environ.get("KINDS", "ins,del").split(",")
FORCE = os.environ.get("FORCE", "0") == "1"
NPZ = f"{BL.OUT}/sv_burden_pool_null{BL.TAG}.npz"
CSV = f"{BL.OUT}/sv_burden_pool_null{BL.TAG}.csv"


def main():
    D = pd.read_csv(f"{BL.OUT}/founder_burden{BL.TAG}.csv", dtype={"founder": str})
    K = np.load(f"{BL.OUT}/founder_kinship.npy")
    Zf = np.load(BL.EF, allow_pickle=True)
    assert np.array_equal(Zf["founders"].astype(str), D.founder.to_numpy())
    Y = np.log((Zf["H0"][None, :] + Zf["DH_rel"] + BL.EPS) / (Zf["H0"][None, :] + BL.EPS)); x = Zf["bio1"]
    L = BL.LMM(K); z = BL.z
    cov = np.column_stack([np.ones(len(D)), z(D.callrate_sv), z(D.callrate_snp), D.is_cactus.astype(float)])
    rate = lambda c: (D[c] / D[f"called_{c}"]).to_numpy(float)

    def fit_all(t, ref):
        X = np.column_stack([cov, z(t), z(ref)]); b = np.empty(Y.shape[0]); v = np.empty(Y.shape[0])
        for s in range(Y.shape[0]):
            bb, se, _ = L.fit(Y[s], X); b[s] = bb[-2]; v[s] = se[-2] ** 2
        w = 1 / v; pooled = (w * b).sum() / w.sum(); pse = np.sqrt(1 / w.sum())
        Xm = np.column_stack([np.ones(b.size), x]); XtW = Xm.T * w; slope = np.linalg.solve(XtW @ Xm, XtW @ b)[1]
        bm, sem, _ = L.fit(Y.mean(0), X)
        return b, np.sqrt(v), pooled, pse, slope, bm[-2], sem[-2]

    old = dict(np.load(NPZ)) if os.path.exists(NPZ) else {}
    out = dict(old); rows = []
    for kind in KINDS:
        if not FORCE and f"{kind}_null_pooled" in old and old[f"{kind}_null_pooled"].size >= NPERM:
            print(f"[{kind}] already in {NPZ} ({old[f'{kind}_null_pooled'].size} draws) -- skipped (FORCE=1 to redo)")
        else:
            t = rate(f"{kind}_shared_CDS_all"); ref = rate(f"{kind}_shared_intergenic_all")
            b, se, pooled, pse, slope, bm, sem = fit_all(t, ref)
            fit = np.polyval(np.polyfit(ref, t, 1), ref); res = t - fit
            rng = np.random.default_rng(0); t0 = time.time()
            NB = np.empty((NPERM, Y.shape[0])); NPO = np.empty(NPERM); NS = np.empty(NPERM); NM = np.empty(NPERM)
            for i in range(NPERM):
                bb, _, po, _, sl, m_, _ = fit_all(fit + rng.permutation(res), ref)
                NB[i] = bb; NPO[i] = po; NS[i] = sl; NM[i] = m_
            out.update({f"{kind}_beta": b, f"{kind}_se": se, f"{kind}_obs": np.array([pooled, pse, slope, bm, sem]),
                        f"{kind}_null_beta": NB.astype(np.float32), f"{kind}_null_pooled": NPO,
                        f"{kind}_null_slope": NS, f"{kind}_null_mean": NM})
            out["sites"] = Zf["sites"]; out["bio1"] = x
            np.savez_compressed(NPZ, **out)                       # saved after every kind
            print(f"[{kind}] {NPERM} permutations in {time.time() - t0:.0f}s -> {NPZ}", flush=True)
        pooled, pse, slope, bm, sem = out[f"{kind}_obs"]
        NPO, NS, NM, NB = out[f"{kind}_null_pooled"], out[f"{kind}_null_slope"], out[f"{kind}_null_mean"], out[f"{kind}_null_beta"]
        iu = np.triu_indices(NB.shape[1], 1)
        rows.append(dict(kind=kind, n_perm=NPO.size, pooled_beta=pooled, se_independent=pse, se_permutation=NPO.std(),
                         se_ratio=NPO.std() / pse, p_independent=2 * stats.norm.sf(abs(pooled / pse)),
                         p_permutation=(np.abs(NPO) >= abs(pooled)).mean(),
                         mean_fitness_beta=bm, mean_fitness_se=sem, mean_fitness_p=2 * stats.norm.sf(abs(bm / sem)),
                         mean_fitness_p_permutation=(np.abs(NM) >= abs(bm)).mean(),
                         slope_bio1=slope, slope_sd_null=NS.std(), slope_p_permutation=(np.abs(NS) >= abs(slope)).mean(),
                         null_corr_between_gardens=float(np.median(np.corrcoef(NB.T)[iu]))))
    R = pd.DataFrame(rows); R.to_csv(CSV, index=False)
    print(R.round(5).T.to_string()); print(f"[wrote] {CSV}")


if __name__ == "__main__":
    main()
