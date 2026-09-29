#!/usr/bin/env python
"""The full founder burden model: SNPs, insertions and deletions fitted TOGETHER.

  y ~ 1 + snp_CDS + snp_intergenic + ins_CDS + ins_intergenic + del_CDS + del_intergenic + cactus + u,
      u ~ N(0, s2g K)
y = founder log-fitness (mean over the 31 gardens, and per garden), as `_sv_burden_lmm.py`. Burden terms
are the shared-variant rates (carried / called, >= 2 carriers, called in >= 90% of founders), standardised
for the fit. `cactus` = assembly route (the one technical term kept; the two call rates were redundant with
it, see `_sv_burden_model_simplify.py`). Same LMM class and kinship matrix as `_sv_burden_lmm.py`.

Two scales are written for every burden term:
  per SD       the fitted coefficient (one SD of founder-to-founder variation in that rate)
  per variant  coefficient / (SD of the rate x median number of that class's variants called per founder):
               the change in log-fitness for one more carried variant of that class, the others held fixed.

Permutation null (founder labels, Kennedy / residualised-predictor scheme, as the headline permutation of
`_sv_burden_lmm.py`): for each tested CDS term, its rate is regressed on all other terms, the residual is
permuted across founders once per draw, and the same permuted vector is used in the mean-fitness fit and in
all 31 per-garden fits. Gives valid p-values for the mean-fitness coefficient, the pooled per-garden
coefficient and the per-garden slope on garden bio1, keeping the dependence between gardens (same founders).

CONTEXT=promoter|UTR|intron|TE (base model): that context replaces CDS as the tested term; outputs tagged `_{CONTEXT}`.
CONTEXT=CDS,promoter: both contexts tested together, each class = coding + promoter + intergenic; tag `_CDS+promoter`.
MODEL=synsplit: SNPs entered as nonsynonymous + synonymous coding burden (`_sv_burden_snp_syn.py`) instead of
coding + intergenic; outputs tagged `_synsplit`, the base model's files are untouched.
Env: kmate, SV_UNIT=event. NPERM (default 3000). Cached: a term already in the npz with >= NPERM draws is
not redone unless FORCE=1.
Writes results/sv_adaptive/sv_burden_full_model{TAG}.csv          (mean fitness, all terms, both scales)
       results/sv_adaptive/sv_burden_full_model_by_garden{TAG}.csv (per garden, all terms, both scales)
       results/sv_adaptive/sv_burden_full_model_perm{TAG}.npz      (null draws per tested term)
"""
import os, sys, time
import numpy as np, pandas as pd
from scipy import stats

os.environ.setdefault("SV_UNIT", "event")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _sv_burden_lmm as BL

NPERM = int(os.environ.get("NPERM", 3000))
FORCE = os.environ.get("FORCE", "0") == "1"
MODEL = os.environ.get("MODEL", "base")          # base | synsplit
if MODEL == "synsplit":
    # SNPs as nonsynonymous + synonymous coding burden (`_sv_burden_snp_syn.py`): synonymous is the SNP class's
    # within-gene neutral baseline and the negative control; SVs keep coding + intergenic
    TERMS = ["snpnon_CDS", "snpsyn_CDS", "ins_CDS", "ins_intergenic", "del_CDS", "del_intergenic"]
    TESTED = ["snpnon_CDS", "snpsyn_CDS", "ins_CDS", "del_CDS"]
else:
    CONTEXT = os.environ.get("CONTEXT", "CDS")      # tested context(s), comma-separated; each class also gets its intergenic rate
    CTXS = CONTEXT.split(",")
    KINDS = ("snp", "ins", "del")
    TERMS = [f"{k}_{x}" for k in KINDS for x in CTXS + ["intergenic"]]
    TESTED = [f"{k}_{x}" for k in KINDS for x in CTXS]
OUT = BL.OUT; TAG = BL.TAG + ("_synsplit" if MODEL == "synsplit" else "" if CONTEXT == "CDS" else "_" + "+".join(CONTEXT.split(",")))
F_MEAN = f"{OUT}/sv_burden_full_model{TAG}.csv"
F_GARD = f"{OUT}/sv_burden_full_model_by_garden{TAG}.csv"
F_PERM = f"{OUT}/sv_burden_full_model_perm{TAG}.npz"


def main():
    D = pd.read_csv(f"{OUT}/founder_burden{BL.TAG}.csv", dtype={"founder": str})
    if MODEL == "synsplit":
        S = pd.read_csv(f"{OUT}/founder_burden_snpsyn.csv", dtype={"founder": str})
        assert np.array_equal(S.founder.to_numpy(), D.founder.to_numpy()); D = pd.concat([D, S.drop(columns="founder")], axis=1)
    K = np.load(f"{OUT}/founder_kinship.npy")
    Zf = np.load(BL.EF, allow_pickle=True)
    assert np.array_equal(Zf["founders"].astype(str), D.founder.to_numpy())
    Y = np.log((Zf["H0"][None, :] + Zf["DH_rel"] + BL.EPS) / (Zf["H0"][None, :] + BL.EPS)); y = Y.mean(0)
    sites, bio1, zone = Zf["sites"].astype(int), Zf["bio1"].astype(float), Zf["zone"].astype(str)
    L = BL.LMM(K); n = len(D)
    raw = {t: (D[f"{t.split('_')[0]}_shared_{t.split('_')[1]}_all"] /
               D[f"called_{t.split('_')[0]}_shared_{t.split('_')[1]}_all"]).to_numpy(float) for t in TERMS}
    nmed = {t: float(np.median(D[f"called_{t.split('_')[0]}_shared_{t.split('_')[1]}_all"])) for t in TERMS}
    sd = {t: raw[t].std() for t in TERMS}
    Z = {t: (raw[t] - raw[t].mean()) / sd[t] for t in TERMS}
    cac = D.is_cactus.astype(float).to_numpy()
    names = TERMS + ["cactus"]

    def design(Zd):
        return np.column_stack([np.ones(n)] + [Zd[t] for t in TERMS] + [cac])

    def fit_set(X):
        b, se, h2 = L.fit(y, X)
        G = [L.fit(Y[s], X) for s in range(Y.shape[0])]
        return b, se, h2, np.array([g[0] for g in G]), np.array([g[1] for g in G])

    def pooled_slope(b, se):
        w = 1 / se ** 2; po = (w * b).sum() / w.sum()
        Xm = np.column_stack([np.ones(b.size), bio1]); XtW = Xm.T * w
        return po, np.sqrt(1 / w.sum()), np.linalg.solve(XtW @ Xm, XtW @ b)[1]

    X = design(Z); b, se, h2, GB, GS = fit_set(X)
    scale = {t: 1 / (sd[t] * nmed[t]) for t in TERMS}; scale["cactus"] = 1.0

    # ---- permutation null for each tested CDS term (cached)
    P = dict(np.load(F_PERM)) if os.path.exists(F_PERM) else {}
    for t in TESTED:
        if not FORCE and f"{t}_null_mean" in P and P[f"{t}_null_mean"].size >= NPERM:
            print(f"[{t}] permutation null cached ({P[f'{t}_null_mean'].size} draws)"); continue
        j = names.index(t) + 1
        Xo = np.delete(X, j, axis=1); fit = Xo @ np.linalg.lstsq(Xo, X[:, j], rcond=None)[0]; res = X[:, j] - fit
        rng = np.random.default_rng(0); t0 = time.time()
        NM = np.empty(NPERM); NG = np.empty((NPERM, Y.shape[0]), np.float32); NPO = np.empty(NPERM); NS = np.empty(NPERM)
        for i in range(NPERM):
            Xp = X.copy(); Xp[:, j] = fit + rng.permutation(res)
            bp, _, _, gb, gs = fit_set(Xp)
            NM[i] = bp[j]; NG[i] = gb[:, j]; po, _, sl = pooled_slope(gb[:, j], gs[:, j]); NPO[i] = po; NS[i] = sl
        P.update({f"{t}_null_mean": NM, f"{t}_null_garden": NG, f"{t}_null_pooled": NPO, f"{t}_null_slope": NS})
        np.savez_compressed(F_PERM, **P)
        print(f"[{t}] {NPERM} permutations in {time.time() - t0:.0f}s -> {F_PERM}", flush=True)

    # ---- mean-fitness table
    rows = []
    for i, t in enumerate(names, 1):
        r = dict(term=t, kind=t.split("_")[0] if t != "cactus" else "technical",
                 context=t.split("_")[1] if t != "cactus" else "assembly route",
                 beta_sd=b[i], se_sd=se[i], p_wald=2 * stats.norm.sf(abs(b[i] / se[i])),
                 beta_per_variant=b[i] * scale[t], se_per_variant=se[i] * scale[t],
                 variants_per_sd=(sd[t] * nmed[t]) if t != "cactus" else np.nan,
                 called_median=nmed.get(t, np.nan), h2=h2)
        po, pse, sl = pooled_slope(GB[:, i], GS[:, i])
        r.update(pooled_sd=po, pooled_se_independent=pse, slope_bio1_sd=sl, negative_in=int((GB[:, i] < 0).sum()))
        if t in TESTED:
            NM, NPO, NS = P[f"{t}_null_mean"], P[f"{t}_null_pooled"], P[f"{t}_null_slope"]
            r.update(p_perm=(np.abs(NM) >= abs(b[i])).mean(), se_perm_mean=NM.std(),
                     pooled_se_perm=NPO.std(), p_perm_pooled=(np.abs(NPO) >= abs(po)).mean(),
                     slope_sd_null=NS.std(), p_perm_slope=(np.abs(NS) >= abs(sl)).mean(), n_perm=NM.size)
        rows.append(r)
    pd.DataFrame(rows).to_csv(F_MEAN, index=False)

    # ---- per-garden table
    gr = []
    for i, t in enumerate(names, 1):
        for s in range(Y.shape[0]):
            gr.append(dict(term=t, site=sites[s], bio1=bio1[s], zone=zone[s], beta_sd=GB[s, i], se_sd=GS[s, i],
                           beta_per_variant=GB[s, i] * scale[t], se_per_variant=GS[s, i] * scale[t],
                           p_wald=2 * stats.norm.sf(abs(GB[s, i] / GS[s, i]))))
    pd.DataFrame(gr).to_csv(F_GARD, index=False)
    T = pd.read_csv(F_MEAN)
    print(T[["term", "beta_sd", "se_sd", "p_wald", "beta_per_variant", "se_per_variant"] +
            [c for c in ("p_perm", "p_perm_pooled", "p_perm_slope") if c in T]].round(6).to_string(index=False))
    print(f"[wrote] {F_MEAN}, {F_GARD}")


if __name__ == "__main__":
    main()
