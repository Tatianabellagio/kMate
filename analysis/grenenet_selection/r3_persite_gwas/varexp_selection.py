#!/usr/bin/env python
"""Variance of the ecotype-selection trait explained by SNP vs non-SNP variation.

The phase-1-vs-kMate question: how much of the per-founder selection response (logit-slope s of
founder frequency over gens 0..3, per site) is captured by SNPs vs by non-SNP variation (indels +
SVs), and how much do we GAIN by adding the non-SNP layer? Three complementary estimators, each
on the analyzable-founder subset (present at founding), traits = per-site s + climate-zone/global
aggregate axes:

  (1) MARGINAL single-GRM REML h2  -- how much each class explains ALONE (K_snp, K_nonsnp, K_indel,
      K_sv, K_all, K_snp_matched). Because non-SNPs are in LD with SNPs these OVERLAP
      (h2_snp + h2_nonsnp > h2_all); the overlap is itself a result.
  (2) JOINT 2-GRM REML (K_snp + K_nonsnp) -- partition into SNP / non-SNP / shared variance, with a
      likelihood-ratio test for whether K_nonsnp adds anything beyond K_snp. Ill-conditioned (the two
      GRMs are correlated); report with SE.
  (3) GBLUP cross-validated prediction -- leave-fold-out predictive R^2 for SNP-only, non-SNP-only,
      both; GAIN = R2(both) - R2(snp). The robust "how much do we gain" number at N=212.

Fairness: K_snp_matched (SNPs subsampled to the non-SNP count & MAC spectrum) isolates whether any
SNP edge is just "more markers". Reports marginal (each class alone) AND conditional gain.

Reads analysis/grenenet_selection/varexp/{selection_s_matrix.npz, class_grms.npz}. Writes varexp.csv,
varexp_persite.npz, varexp_meta.json. Env: kmate. Light (all modeling on 212x212 from cache).
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/varexp"
CLASSES = ["snp", "nonsnp", "indel", "sv", "all", "snp_matched"]
NFOLD, NREP, SEED = 6, 10, 0


# ---------------------------------------------------------------- REML (AI) --
def aireml(y, X, Ks, n_iter=100, tol=1e-4):
    """Average-Information REML for V = sum_k s_k^2 K_k + s_e^2 I (I is the last 'GRM').

    Returns (theta[len(Ks)+1] variance comps, cov[q x q], loglik). Falls back to an EM step when
    an AI update leaves the boundary. Small n (<=231) so dense solves are cheap."""
    n = len(y); comps = Ks + [np.eye(n)]; q = len(comps)
    vy = np.var(y)
    theta = np.full(q, vy / q)
    ll_prev = -np.inf
    for _ in range(n_iter):
        V = sum(t * C for t, C in zip(theta, comps))
        Vi = np.linalg.inv(V + 1e-8 * np.eye(n))
        XtViX = X.T @ Vi @ X
        XtViX_i = np.linalg.inv(XtViX)
        P = Vi - Vi @ X @ XtViX_i @ X.T @ Vi
        Py = P @ y
        # loglik (REML)
        s_ld = np.linalg.slogdet(V + 1e-8 * np.eye(n))[1]
        x_ld = np.linalg.slogdet(XtViX)[1]
        ll = -0.5 * (s_ld + x_ld + y @ Py)
        PC = [P @ C for C in comps]
        score = np.array([-0.5 * (np.trace(PC[i]) - Py @ comps[i] @ Py) for i in range(q)])
        AI = np.empty((q, q))
        for i in range(q):
            for j in range(i, q):
                AI[i, j] = AI[j, i] = 0.5 * (Py @ comps[i] @ P @ comps[j] @ Py)
        try:
            step = np.linalg.solve(AI + 1e-8 * np.eye(q), score)
        except np.linalg.LinAlgError:
            step = np.zeros(q)
        new = theta + step
        if np.any(new < 0):                      # EM fallback keeps components >= 0
            new = np.array([max(theta[i] + theta[i] ** 2 * score[i] * 2 / n, 1e-9)
                            for i in range(q)])
        theta = np.maximum(new, 1e-9)
        if abs(ll - ll_prev) < tol:
            break
        ll_prev = ll
    try:
        cov = np.linalg.inv(AI + 1e-8 * np.eye(q))
    except np.linalg.LinAlgError:
        cov = np.full((q, q), np.nan)
    return theta, cov, ll


def h2_1grm(y, K):
    """Marginal h2 for one GRM + delta-method SE + LRT p vs h2=0."""
    n = len(y); X = np.ones((n, 1))
    th, cov, ll = aireml(y, X, [K])
    vg, ve = th
    mbar = np.mean(np.diag(K))
    tot = vg * mbar + ve
    h2 = vg * mbar / tot
    # delta-method SE of h2 = (vg*mbar)/(vg*mbar+ve)
    g = np.array([ve * mbar, -vg * mbar]) / tot ** 2      # d h2 / d(vg,ve)
    se = float(np.sqrt(max(g @ cov @ g, 0)))
    # LRT vs null (no GRM): fit intercept-only variance
    ll0 = -0.5 * (n * np.log(2 * np.pi * np.var(y)) + n)
    lr = max(2 * (ll - ll0), 0)
    p = 0.5 * stats.chi2.sf(lr, 1)                          # 1-df boundary mixture
    return h2, se, p


def joint_2grm(y, K1, K2):
    """Joint partition of variance into K1 / K2 / residual + LRT for adding K2 over K1."""
    n = len(y); X = np.ones((n, 1))
    th, cov, ll2 = aireml(y, X, [K1, K2])
    _, _, ll1 = aireml(y, X, [K1])
    m1, m2 = np.mean(np.diag(K1)), np.mean(np.diag(K2))
    scaled = np.array([th[0] * m1, th[1] * m2, th[2]])
    frac = scaled / scaled.sum()
    lr = max(2 * (ll2 - ll1), 0)
    p_add = 0.5 * stats.chi2.sf(lr, 1)
    return frac, p_add           # frac = [K1, K2, resid]


# ---------------------------------------------------- GBLUP cross-validation --
def gblup_cv(y, Ks, nfold=NFOLD, nrep=NREP, seed=SEED):
    """Repeated k-fold predictive correlation & R^2 for the combined GRM list Ks."""
    n = len(y); X = np.ones((n, 1))
    preds = np.zeros(n); cnt = np.zeros(n)
    rng = np.random.default_rng(seed)
    for rep in range(nrep):
        idx = rng.permutation(n); folds = np.array_split(idx, nfold)
        for te in folds:
            tr = np.setdiff1d(np.arange(n), te)
            Ks_tr = [K[np.ix_(tr, tr)] for K in Ks]
            th, _, _ = aireml(y[tr], X[tr], Ks_tr)
            V = sum(th[k] * Ks[k][np.ix_(tr, tr)] for k in range(len(Ks))) + th[-1] * np.eye(len(tr))
            Vi = np.linalg.inv(V + 1e-8 * np.eye(len(tr)))
            xo = np.ones(len(tr))
            beta = float((xo @ Vi @ y[tr]) / (xo @ Vi @ xo))
            resid = y[tr] - beta
            Kcross = sum(th[k] * Ks[k][np.ix_(te, tr)] for k in range(len(Ks)))
            preds[te] += beta + Kcross @ Vi @ resid; cnt[te] += 1
    yhat = preds / np.maximum(cnt, 1)
    r = float(np.corrcoef(y, yhat)[0, 1])
    r2 = 1.0 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)
    return r, float(r2)


def rint(x):
    r = stats.rankdata(x); return stats.norm.ppf((r - 0.5) / len(x))


def main():
    S = np.load(f"{OUT}/selection_s_matrix.npz", allow_pickle=True)
    G = np.load(f"{OUT}/class_grms.npz", allow_pickle=True)
    nmark = json.loads(str(G["n_markers"]))
    tf = S["founders"].astype(str); gf = G["founders"].astype(str)
    ana = S["analyzable"].astype(bool)
    # align GRM founders -> trait founders, keep analyzable & present in both
    gpos = {f: i for i, f in enumerate(gf)}
    keep = np.array([ana[i] and tf[i] in gpos for i in range(len(tf))])
    order = [gpos[tf[i]] for i in range(len(tf)) if keep[i]]
    Ks = {c: G[f"K_{c}"][np.ix_(order, order)] for c in CLASSES}
    n = len(order)
    Smat = S["S"][:, keep]                       # n_site x n
    sites = S["sites"]; bio1 = S["bio1"]
    print(f"founders analyzable & genotyped: {n}   sites: {len(sites)}", flush=True)
    print("n_markers:", nmark, flush=True)
    print("corr(K_snp,K_nonsnp) offdiag = %.3f" %
          np.corrcoef(Ks["snp"][np.triu_indices(n, 1)],
                      Ks["nonsnp"][np.triu_indices(n, 1)])[0, 1], flush=True)

    # ---- aggregate axes: global mean + bio1-tercile zone means ----
    zone = pd.qcut(bio1, 3, labels=["cold", "mid", "hot"]).astype(str)
    axes = {"w_global": Smat.mean(0)}
    for z in ("cold", "mid", "hot"):
        axes[f"w_{z}"] = Smat[zone == z].mean(0)

    rows = []

    def run_trait(name, ttype, y_raw):
        for transform, y in (("raw", y_raw), ("rint", rint(y_raw))):
            y = np.asarray(y, float)
            for c in CLASSES:                                     # (1) marginal
                h2, se, p = h2_1grm(y, Ks[c])
                rows.append(dict(trait=name, ttype=ttype, transform=transform,
                                 method="marginal", cls=c, value=h2, se=se, p=p))
            frac, p_add = joint_2grm(y, Ks["snp"], Ks["nonsnp"])  # (2) joint
            for cl, v in zip(["snp", "nonsnp", "resid"], frac):
                rows.append(dict(trait=name, ttype=ttype, transform=transform,
                                 method="joint2", cls=cl, value=v, se=np.nan,
                                 p=(p_add if cl == "nonsnp" else np.nan)))
            r_s, r2_s = gblup_cv(y, [Ks["snp"]])                  # (3) prediction CV
            r_n, r2_n = gblup_cv(y, [Ks["nonsnp"]])
            r_b, r2_b = gblup_cv(y, [Ks["snp"], Ks["nonsnp"]])
            r_m, r2_m = gblup_cv(y, [Ks["snp_matched"]])
            for cl, rr, r2 in [("snp", r_s, r2_s), ("nonsnp", r_n, r2_n),
                               ("both", r_b, r2_b), ("snp_matched", r_m, r2_m)]:
                rows.append(dict(trait=name, ttype=ttype, transform=transform,
                                 method="predcv", cls=cl, value=r2, se=np.nan, p=rr))  # p col=pearson r
            rows.append(dict(trait=name, ttype=ttype, transform=transform,
                             method="predcv", cls="gain_nonsnp", value=r2_b - r2_s,
                             se=np.nan, p=r_b - r_s))

    for name, y in axes.items():
        print(f"axis {name} ...", flush=True); run_trait(name, "axis", y)

    # ---- per-site scan (raw only, marginal + predcv gain) for the across-site distribution ----
    persite = {c: [] for c in ["h2_snp", "h2_nonsnp", "h2_all",
                               "r2_snp", "r2_nonsnp", "r2_both", "gain"]}
    for i, site in enumerate(sites):
        y = Smat[i]
        for c, key in [("snp", "h2_snp"), ("nonsnp", "h2_nonsnp"), ("all", "h2_all")]:
            persite[key].append(h2_1grm(y, Ks[c])[0])
        _, r2s = gblup_cv(y, [Ks["snp"]]); _, r2n = gblup_cv(y, [Ks["nonsnp"]])
        _, r2b = gblup_cv(y, [Ks["snp"], Ks["nonsnp"]])
        persite["r2_snp"].append(r2s); persite["r2_nonsnp"].append(r2n)
        persite["r2_both"].append(r2b); persite["gain"].append(r2b - r2s)
        rows.append(dict(trait=f"site{site}", ttype="site", transform="raw",
                         method="predcv", cls="gain_nonsnp", value=r2b - r2s, se=np.nan, p=np.nan))
        print(f"  site {site}: h2_snp={persite['h2_snp'][-1]:.2f} "
              f"h2_nonsnp={persite['h2_nonsnp'][-1]:.2f} gain_r2={r2b-r2s:+.3f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/varexp.csv", index=False)
    np.savez(f"{OUT}/varexp_persite.npz", sites=sites, bio1=bio1,
             **{k: np.array(v) for k, v in persite.items()})

    # ---- summary ----
    ax = df[(df.ttype == "axis") & (df.transform == "raw")]
    def cell(m, c, tr="w_global"):
        v = ax[(ax.method == m) & (ax.cls == c) & (ax.trait == tr)].value
        return float(v.iloc[0]) if len(v) else np.nan
    ps = {k: np.array(v) for k, v in persite.items()}
    summary = {
        "n_founders": n, "n_sites": int(len(sites)), "n_markers": nmark,
        "corr_Ksnp_Knonsnp": float(np.corrcoef(Ks["snp"][np.triu_indices(n, 1)],
                                                Ks["nonsnp"][np.triu_indices(n, 1)])[0, 1]),
        "w_global_marginal_h2": {c: cell("marginal", c) for c in CLASSES},
        "w_global_joint": {c: cell("joint2", c) for c in ["snp", "nonsnp", "resid"]},
        "w_global_predR2": {c: cell("predcv", c) for c in ["snp", "nonsnp", "both",
                                                           "snp_matched", "gain_nonsnp"]},
        "persite_median": {k: round(float(np.median(v)), 3) for k, v in ps.items()},
        "persite_gain_frac_positive": float(np.mean(ps["gain"] > 0)),
    }
    json.dump(summary, open(f"{OUT}/varexp_meta.json", "w"), indent=2)
    print("\n===== SUMMARY (w_global, raw) =====")
    print("marginal h2:", {c: round(cell('marginal', c), 3) for c in CLASSES})
    print("joint2 frac:", {c: round(cell('joint2', c), 3) for c in ['snp', 'nonsnp', 'resid']})
    print("pred R2   :", {c: round(cell('predcv', c), 3)
                          for c in ['snp', 'nonsnp', 'both', 'snp_matched', 'gain_nonsnp']})
    print(f"per-site median gain R2 = {np.median(ps['gain']):+.3f}  "
          f"(frac sites gain>0 = {np.mean(ps['gain'] > 0):.2f})")
    print(f"wrote {OUT}/varexp.csv + varexp_persite.npz + varexp_meta.json")


if __name__ == "__main__":
    main()
