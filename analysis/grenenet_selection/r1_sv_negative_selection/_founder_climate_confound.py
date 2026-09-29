#!/usr/bin/env python
"""Is the climate-graded SV purging signal SV load, or founder local adaptation?

WHY THIS EXISTS. In global mode per-variant AF is AF_v = h.V_pa[:,v] / h.V_called[:,v] --
a deterministic projection of the per-sample founder mixture h. So the per-variant climate
slope beta_v (_compute_s_climate_slope.py) carries information ONLY through (i) which
founders carry v and (ii) how those founders respond to climate. `_founder_load_test.py`
established the founder-level fact that follows from this: founders carrying a higher
FRACTION of insertions decline more in hot/arid gardens (T1).

That reframes the question but does not answer it. "Insertion-rich founders are purged in
hot gardens" has (at least) four readings, and they are NOT distinguished by anything in
the thread so far:

  (H1) LOAD        insertion load is deleterious and more so under heat/drought stress.
  (H2) PROVENANCE  insertion-rich founders happen to come from cold/wet homes, so they
                   lose in hot/arid gardens by ordinary local adaptation. Nothing to do
                   with SVs -- insertion load is just a marker of where a founder is from.
  (H3) PANEL       insertion load tracks cactus-vs-PanGenie panel membership (only 80
                   founders have long-read assemblies), and the EM has a known +41%
                   cactus h-bias (BACKGROUND.md limitation 2). Then "insertion load" is
                   partly an assay covariate, not a genome property.
  (H4) STRUCTURE   231 founders are ~10-25 clades, so a Spearman over 231 points has far
                   fewer effective df than it looks. The correlation may not be
                   significant once kinship is accounted for.

This script tests H2/H3/H4 head-on against H1, at the FOUNDER level, where the causal
question actually lives. It deliberately uses no per-variant temporal data.

Outputs: results/sv_adaptive/founder_climate_confound.{npz,csv} + printed report.
Env: kmate. Reads selection_s_matrix.npz, founder_load_test.npz, class_grms.npz,
the 1001g ecotype bioclim table, and data/founder_split_cactus_pg.json.
"""
import os, sys, json, glob
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

PROJ = lib.PROJ
GEA = lib.GEA
PM = f"{GEA}/common/results/pool_matrices"
OUT = f"{GEA}/r1_sv_negative_selection/results/sv_adaptive"
ECO = ("/global/scratch/users/tbellg/gea_grene-net/key_files/"
       "1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv")
CVARS = ("bio1", "bio18")


# ---------------------------------------------------------------- mixed model
def emma_reml(y, X, K, n_grid=100):
    """Single-variance-component REML: y ~ N(Xb, sg2 * (K + delta I)).

    Returns (beta, se, delta, sg2, loglik). Eigendecomposition-based, n=231 so exact.
    """
    n, p = X.shape
    d, U = np.linalg.eigh(K)
    d = np.maximum(d, 0.0)
    yt = U.T @ y
    Xt = U.T @ X

    def fit(delta):
        h = d + delta
        Xh = Xt / h[:, None]
        A = Xt.T @ Xh
        Ainv = np.linalg.pinv(A)
        beta = Ainv @ (Xh.T @ yt)
        r = yt - Xt @ beta
        rss = float((r * r / h).sum())
        sg2 = rss / (n - p)
        sign, logdetA = np.linalg.slogdet(A)
        ll = -0.5 * ((n - p) * np.log(sg2) + np.log(h).sum() + logdetA + (n - p))
        return beta, sg2, ll, Ainv

    grid = np.logspace(-5, 5, n_grid)
    best = max(((fit(dl), dl) for dl in grid), key=lambda t: t[0][2])
    (beta, sg2, ll, Ainv), delta = best
    se = np.sqrt(np.maximum(np.diag(Ainv) * sg2, 0.0))
    return beta, se, delta, sg2, ll


def partial_spearman(x, y, Z):
    """Spearman of x,y after linearly residualizing the RANKS of both on ranks of Z."""
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    RZ = np.column_stack([stats.rankdata(Z[:, j]) for j in range(Z.shape[1])])
    A = np.column_stack([np.ones(len(rx)), RZ])
    bx, *_ = np.linalg.lstsq(A, rx, rcond=None)
    by, *_ = np.linalg.lstsq(A, ry, rcond=None)
    ex = rx - A @ bx; ey = ry - A @ by
    r = float(np.corrcoef(ex, ey)[0, 1])
    dof = len(rx) - 2 - Z.shape[1]
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    return r, float(2 * stats.t.sf(abs(t), dof))


def zs(v):
    v = np.asarray(v, float)
    return (v - np.nanmean(v)) / np.nanstd(v)


def main():
    # ------------------------------------------------ founder climate response gamma_f
    Z = np.load(f"{GEA}/r3_persite_gwas/results/varexp/selection_s_matrix.npz", allow_pickle=True)
    S = Z["S"]                                  # [sites x founders] founder selection coef
    founders = Z["founders"].astype("U6")
    sites = Z["sites"].astype(int)
    p0f = Z["p0"].astype(np.float64)
    analyz = Z["analyzable"].astype(bool)

    cmeta = pd.concat([pd.read_csv(f"{PM}/pool_gen{g}_snp.meta.csv") for g in (1, 2, 3)],
                      ignore_index=True).groupby("site")[list(CVARS)].mean()
    clim = {k: np.array([cmeta.loc[s, k] if s in cmeta.index else np.nan for s in sites])
            for k in CVARS}

    gamma = {}
    for k in CVARS:
        c = clim[k]; ok = np.isfinite(c)
        g = np.full(founders.size, np.nan)
        for f in range(founders.size):
            col = S[ok, f]; m = np.isfinite(col)
            if m.sum() < 5:
                continue
            cx = c[ok][m] - c[ok][m].mean()
            g[f] = float((cx @ col[m]) / (cx @ cx))
        gamma[k] = g
        print(f"[gamma] {k}: n={np.isfinite(g).sum()} median={np.nanmedian(g):+.5f} sd={np.nanstd(g):.5f}")

    # ------------------------------------------------ founder SV load (from T1 artifact)
    F = np.load(f"{OUT}/founder_load_test.npz", allow_pickle=True)
    assert list(F["founders"].astype("U6")) == list(founders), "founder order mismatch"
    ins_load, del_load, tot = F["ins_load"], F["del_load"], F["tot_carry"]
    ins_frac = ins_load / np.maximum(tot, 1)
    del_frac = del_load / np.maximum(tot, 1)

    # ------------------------------------------------ founder HOME climate (H2)
    eco = pd.read_csv(ECO)
    eco["ecotypeid"] = eco["ecotypeid"].astype(str)
    emap = eco.set_index("ecotypeid")
    home = {}
    for k in CVARS:
        home[k] = np.array([float(emap.loc[f, k]) if f in emap.index else np.nan
                            for f in founders])
    n_home = np.isfinite(home["bio1"]).sum()
    print(f"[home]  matched {n_home}/{founders.size} founders to 1001g bioclim")

    # ------------------------------------------------ panel membership (H3)
    split = json.load(open(f"{PROJ}/data/founder_split_cactus_pg.json"))
    cac = set(str(x) for x in split["cactus"])
    is_cactus = np.array([f in cac for f in founders], dtype=float)
    print(f"[panel] cactus={int(is_cactus.sum())} PG={int((1 - is_cactus).sum())}")

    # ------------------------------------------------ kinship (H4)
    G = np.load(f"{GEA}/r3_persite_gwas/results/varexp/class_grms.npz", allow_pickle=True)
    assert list(G["founders"].astype("U6")) == list(founders), "GRM founder order mismatch"
    K = G["K_snp"].astype(np.float64)
    K = K / np.mean(np.diag(K))

    ok = (analyz & np.isfinite(gamma["bio1"]) & np.isfinite(home["bio1"])
          & np.isfinite(home["bio18"]) & (tot > 0))
    print(f"[keep]  {int(ok.sum())} founders with gamma + load + home climate\n")

    # effective df under kinship: (sum d)^2 / sum d^2 on the kept submatrix
    Ks = K[np.ix_(ok, ok)]
    dK = np.linalg.eigvalsh(Ks); dK = np.maximum(dK, 0)
    n_eff = float(dK.sum() ** 2 / (dK ** 2).sum())
    print(f"[H4] kinship effective df ~ {n_eff:.1f} independent founders "
          f"(vs {int(ok.sum())} nominal)\n")

    rows = []
    print("=" * 78)
    print("H2/H3  Is insertion load confounded with founder PROVENANCE or PANEL?")
    for lbl, v in (("home bio1", home["bio1"]), ("home bio18", home["bio18"]),
                   ("is_cactus", is_cactus), ("founder p0", p0f)):
        r, p = stats.spearmanr(ins_frac[ok], v[ok])
        print(f"  Spearman(ins_frac, {lbl:11s}) = {r:+.3f}  p={p:.2e}")
        rows.append(dict(test="ins_frac_vs", target=lbl, rho=r, p=p))
    ic = is_cactus[ok].astype(bool)
    mw = stats.mannwhitneyu(ins_frac[ok][ic], ins_frac[ok][~ic])
    print(f"  ins_frac  cactus median={np.median(ins_frac[ok][ic]):.4f} vs "
          f"PG median={np.median(ins_frac[ok][~ic]):.4f}  MWU p={mw.pvalue:.2e}")

    print("\n" + "=" * 78)
    print("H1 vs H2  Raw and PARTIAL correlations of the founder climate response gamma")
    for k in CVARS:
        g = gamma[k]
        print(f"\n  --- gamma_{k} ---")
        for lbl, v in (("ins_frac", ins_frac), ("del_frac", del_frac),
                       ("home bio1", home["bio1"]), ("home bio18", home["bio18"]),
                       ("is_cactus", is_cactus), ("founder p0", p0f)):
            r, p = stats.spearmanr(g[ok], v[ok])
            print(f"    Spearman(gamma, {lbl:11s})      = {r:+.3f}  p={p:.2e}")
            rows.append(dict(test=f"gamma_{k}_vs", target=lbl, rho=r, p=p))
        # partial: does ins_frac survive controlling for provenance / panel / p0?
        for lbl, Zc in (
            ("| home bio1+bio18", np.column_stack([home["bio1"][ok], home["bio18"][ok]])),
            ("| is_cactus", is_cactus[ok][:, None]),
            ("| p0", p0f[ok][:, None]),
            ("| home+cactus+p0", np.column_stack([home["bio1"][ok], home["bio18"][ok],
                                                  is_cactus[ok], p0f[ok]])),
        ):
            r, p = partial_spearman(ins_frac[ok], g[ok], Zc)
            print(f"    partial(gamma, ins_frac {lbl:18s}) = {r:+.3f}  p={p:.2e}")
            rows.append(dict(test=f"gamma_{k}_partial", target=lbl, rho=r, p=p))

    print("\n" + "=" * 78)
    print("H4  Founder mixed model  gamma ~ covariates + u,  u ~ N(0, sg2 K_snp)")
    mm = []
    for k in CVARS:
        y = gamma[k][ok]
        specs = [
            ("ins_frac only", [("ins_frac", ins_frac)]),
            ("home only", [("home_bio1", home["bio1"]), ("home_bio18", home["bio18"])]),
            ("ins_frac + home", [("ins_frac", ins_frac), ("home_bio1", home["bio1"]),
                                 ("home_bio18", home["bio18"])]),
            ("full", [("ins_frac", ins_frac), ("home_bio1", home["bio1"]),
                      ("home_bio18", home["bio18"]), ("is_cactus", is_cactus),
                      ("p0", p0f)]),
        ]
        for name, cov in specs:
            X = np.column_stack([np.ones(int(ok.sum()))] + [zs(v[ok]) for _, v in cov])
            beta, se, delta, sg2, ll = emma_reml(zs(y), X, Ks)
            print(f"\n  gamma_{k}  [{name}]   (delta={delta:.3g})")
            for j, (nm, _) in enumerate(cov, start=1):
                zst = beta[j] / se[j] if se[j] > 0 else np.nan
                pv = 2 * stats.norm.sf(abs(zst))
                flag = "***" if pv < 1e-3 else "**" if pv < 0.01 else "*" if pv < 0.05 else "ns"
                print(f"      {nm:11s} beta={beta[j]:+.3f} se={se[j]:.3f} "
                      f"z={zst:+.2f} p={pv:.3g} {flag}")
                mm.append(dict(axis=k, model=name, term=nm, beta=beta[j], se=se[j],
                               z=zst, p=pv))

    pd.DataFrame(rows).to_csv(f"{OUT}/founder_climate_confound_corr.csv", index=False)
    pd.DataFrame(mm).to_csv(f"{OUT}/founder_climate_confound_mm.csv", index=False)
    np.savez_compressed(f"{OUT}/founder_climate_confound.npz",
                        founders=founders, keep=ok, ins_frac=ins_frac, del_frac=del_frac,
                        gamma_bio1=gamma["bio1"], gamma_bio18=gamma["bio18"],
                        home_bio1=home["bio1"], home_bio18=home["bio18"],
                        is_cactus=is_cactus, p0f=p0f, n_eff=n_eff)
    print(f"\n[wrote] {OUT}/founder_climate_confound.npz + _corr.csv + _mm.csv")


if __name__ == "__main__":
    main()
