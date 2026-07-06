#!/usr/bin/env python
"""Predict ecotype performance: kinship ("who you are") vs climate-match (local adaptation).

NEW, standalone (does NOT touch any other agent's code). Answers the user's question: can we
predict a founder's selection response at a common garden from WHO IT IS (kinship / generalist
genetic merit) alone, or does CLIMATE MATCH (origin climate x garden climate) add predictive
power -> evidence of climate-driven local adaptation *beyond* clade.

Data
  performance  S[f,site] = per-founder selection = mean-over-plots logit-slope of the founder's
    genome-wide frequency h_f over generations 0..3 at each persistent site (== cross_site_winners'
    Smat). Cached to results/grenenet_gea/hapfreq/cross_site_S_matrix.npz (heavy Lustre I/O; build
    once via sbatch, then modeling is instant).
  kinship  K = founder genomic-relationship matrix from build_genotype (231 x 231).
  origin   o_f = founder home WorldClim bio1 (gea_grene-net worldclim_ecotypesdata; 231/231 matched
    by ecotype id). Multivariate origin climate (bio1..19, standardized) also loaded.
  garden   c_site = site ERA5 bio1 (gea_grene-net bioclimvars_experimental_sites_era5).
  match    m[f,site] = -(z(o_f) - z(c_site))^2   (home-field advantage; gamma>0 = local adaptation)
           and a multivariate version -||z(O_f) - z(C_site)||^2 over all 19 bioclim vars.

Model (GxE reaction norm / variance partition)
  S = mu + site_effect + gamma*match + u_f + e,   u_f ~ N(0, sigma_g^2 K)   [kinship random effect]
  site fixed effects absorb each garden's mean (incl. garden-climate main effect); u_f (constant
  across sites) is the GENERALIST "who you are" value -> ONE ranking everywhere; the match term is
  the ONLY thing that re-ranks founders for a new climate.

Decisive test = LEAVE-ONE-SITE-OUT CV
  A) kinship-only GBLUP : generalist value g_f from the other sites -> same prediction at held-out garden.
  B) kinship + match    : g_f + gamma*match(o_f, c_heldout) -> re-ranks per garden.
  Metric = correlation(predicted, observed founder performance) at the held-out garden, averaged.
  B>A  => climate adds predictive power for a NEW garden = adaptation, not just "who you are".
Plus: in-sample REML variance partition (sigma_g^2 vs match) and a permutation test on gamma
(shuffle founder origin-climate labels -> preserves kinship & marginal performance, breaks match).

Writes results/grenenet_gea/hapfreq/predict_ecotype_performance.{csv,json,png}. Env: kmate.
"""
import os, sys, json, glob
import numpy as np
import pandas as pd
from scipy import stats, optimize
from scipy.linalg import cho_factor, cho_solve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from founder_genotype import build_genotype
from ecotype_selection_site import genome_h

H = "results/grenenet_gea/hapfreq"
GEA = "/global/scratch/users/tbellg/gea_grene-net"
WIN = "results/grenenet_kmate_window"
SEED = "results/grenenet_kmate_window_seedmix"
Tg = np.array([0.0, 1.0, 2.0, 3.0]); EPS = 1e-3
BIOS = [f"bio{i}" for i in range(1, 20)]


def logit(p):
    p = np.clip(p, EPS, 1 - EPS); return np.log(p / (1 - p))


EXCLUDE_SITES = {int(x) for x in os.environ.get("EXCLUDE_SITES", "33").split(",") if x.strip()}


def build_S(founders):
    """per-founder selection S[site,f] = mean-over-plots logit-frequency slope of genome-wide h.

    The strongest selection is in generation 1, so we do NOT require a site to persist to gen 3:
    per plot we fit the OLS slope of logit(h) over whatever timepoints are available, always
    anchored at the founding seed mix (t=0) and requiring at least gen 1 (seedmix->gen1 minimum).
    Sites with the full 0..3 trajectory use all of it; gen-1-only sites use seedmix->gen1. A site
    is kept if it has >=2 such plots. Sites in EXCLUDE_SITES (default 33, an outlier garden) are
    dropped. Cached separately from the strict gen-1-2-3 matrix."""
    cache = f"{H}/cross_site_S_matrix_gen1.npz"
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return z["S"], list(z["sites"]), z["founders"]
    seeds = sorted({p.split("/")[-1].split("_Chr")[0]
                    for p in glob.glob(f"{SEED}/*_Chr1.h_blocks_per_chrom.npz")})
    p0 = np.mean([genome_h(s, SEED) for s in seeds], 0)
    lp0 = logit(p0)
    pt = lib.pool_table()
    pt = pt[pt.sampleid.astype(str).apply(
        lambda x: os.path.exists(f"{WIN}/{x}_Chr1.h_blocks_per_chrom.npz"))]
    gh_cache = {}
    def gh(s):
        if s not in gh_cache:
            gh_cache[s] = genome_h(s, WIN)
        return gh_cache[s]
    S = {}
    for site, sd in pt.groupby("site"):
        if int(site) in EXCLUDE_SITES:
            continue
        cell = {}
        for (gen, plot), g in sd.groupby(["generation", "plot"]):
            hs = [gh(str(x)) for x in g.sampleid]; hs = [h for h in hs if h is not None]
            if hs:
                cell[(int(gen), int(plot))] = np.mean(hs, 0)
        plotgens = {}
        for (gen, plot) in cell:
            plotgens.setdefault(plot, []).append(gen)
        sl = []
        for pl, gens in plotgens.items():
            gens = sorted(g for g in gens if g in (1, 2, 3))       # evolved timepoints present
            if not gens:                                            # need at least gen 1
                continue
            ts = np.array([0.0] + [float(g) for g in gens])         # founding t=0 always anchored
            ys = np.vstack([lp0] + [logit(cell[(g, pl)]) for g in gens])  # (npt, nF)
            tc = ts - ts.mean()
            sl.append((tc[:, None] * ys).sum(0) / (tc ** 2).sum())  # OLS per-generation slope
        if len(sl) < 2:                                             # need >=2 replicate plots
            continue
        S[int(site)] = np.mean(sl, 0)
    sites = sorted(S)
    Smat = np.vstack([S[s] for s in sites])
    np.savez(cache, S=Smat, sites=np.array(sites), founders=np.array(founders))
    print(f"[cache] wrote {cache}  ({Smat.shape[0]} sites x {Smat.shape[1]} founders); "
          f"excluded {sorted(EXCLUDE_SITES)}")
    return Smat, sites, np.array(founders)


def reml_fit(y, X, B, d):
    """REML for y = X beta + Z_f g + e,  g ~ N(0, sg^2 K),  e ~ N(0, se^2 I), the founder effect
    g shared across that founder's sites. Proper low-rank Woodbury (the earlier version wrongly
    assumed the obs-expanded loading was orthonormal; it is not -- each founder appears at ~nS
    sites so B'B ~ nS*I, not I).

    B = Z_f Vk (n x r) is the obs-row-expanded matrix of K's eigenvectors, d = eigvals(K).
    With P = B*sqrt(d):  Z_f K Z_f' = P P',  V_unit = P P' + delta I,  delta = se^2/sg^2.
      V_unit^-1 x   = (x - P (delta I + P'P)^-1 P'x) / delta
      logdet V_unit = (n - r) log delta + logdet(delta I + P'P)
    Profiles sg^2 out; minimizes over log delta. Returns beta, sg^2, se^2, cov(beta), and
    Vir = V_unit^-1 (y - X beta)  (-> founder BLUP  g = K Z_f' Vir,  sg^2 cancels)."""
    n, p = X.shape; r = B.shape[1]
    P = B * np.sqrt(np.clip(d, 0, None))[None, :]      # n x r
    PtP = P.T @ P; PtX = P.T @ X; Pty = P.T @ y
    XtX = X.T @ X; Xty = X.T @ y; yty = float(y @ y)

    def solve(delta):
        M = delta * np.eye(r) + PtP; cf = cho_factor(M, lower=True)
        iMPtX = cho_solve(cf, PtX); iMPty = cho_solve(cf, Pty)
        XtViX = (XtX - PtX.T @ iMPtX) / delta
        XtViy = (Xty - PtX.T @ iMPty) / delta
        yViy = (yty - Pty @ iMPty) / delta
        beta = np.linalg.solve(XtViX, XtViy)
        rVir = yViy - beta @ XtViy                     # = (y-Xb)'Vu^-1(y-Xb) via GLS identity
        ld_M = np.linalg.slogdet(M)[1]
        return beta, XtViX, rVir, ld_M

    def dev(logdelta):
        delta = np.exp(logdelta)
        beta, XtViX, rVir, ld_M = solve(delta)
        logdetVu = (n - r) * np.log(delta) + ld_M
        ldX = np.linalg.slogdet(XtViX)[1]
        return 0.5 * ((n - p) * np.log(max(rVir, 1e-300)) + logdetVu + ldX)

    res = optimize.minimize_scalar(dev, bounds=(-12, 12), method="bounded")
    delta = np.exp(res.x)
    beta, XtViX, rVir, _ = solve(delta)
    sg2 = rVir / (n - p); se2 = delta * sg2
    cov_beta = np.linalg.inv(XtViX) * sg2
    M = delta * np.eye(r) + PtP; cf = cho_factor(M, lower=True)
    rr = y - X @ beta
    Vir = (rr - P @ cho_solve(cf, P.T @ rr)) / delta
    return beta, sg2, se2, cov_beta, Vir


def main():
    G, founders, reg = build_genotype()
    founders = [str(f) for f in founders]
    nF = len(founders)
    S, sites, _ = build_S(founders)               # (nsite, nF)
    nS = len(sites)
    print(f"performance matrix: {nS} sites x {nF} founders")

    # kinship -- from COMMON haplotype clusters only (MAF >= 5%). Without the floor the
    # 1/sqrt(p(1-p)) standardization gives a private haplotype (1 founder, p~0.004) a Z~15 on its
    # lone carrier, so rare/private haplotypes dominate K and inflate pairwise relatedness -- the
    # same failure mode as the founder-GWAS audit (GRM from MAF>=5%, test to lower MAF).
    pall = G.mean(0); common = (pall >= 0.05) & (pall <= 0.95)
    Gc = G[:, common]; pc = Gc.mean(0)
    Z = (Gc - pc) / np.sqrt(pc * (1 - pc) + 1e-9); K = (Z @ Z.T) / Gc.shape[1]
    K /= np.mean(np.diag(K))                       # scale to ~unit diagonal
    print(f"kinship K from {int(common.sum()):,}/{G.shape[1]:,} common haplotype clusters "
          f"(MAF>=5%); dropped {int((~common).sum()):,} rare")

    # origin climate (founder) and garden climate (site)
    wc = pd.read_csv(f"{GEA}/idea_fromind_to_pop/worldclim_ecotypesdata_sorted_20240517.csv", sep="\t")
    wc["ecotypeid"] = wc["ecotypeid"].astype(str)
    wc = wc.set_index("ecotypeid").reindex(founders)
    O = wc[BIOS].to_numpy(float)                   # (nF, 19) origin climate
    o1 = O[:, 0]                                   # origin bio1
    site_clim = pd.read_csv(f"{GEA}/bioclimvars_experimental_sites_era5.csv")
    site_clim["site"] = site_clim["site"].astype(int)
    site_clim = site_clim.set_index("site").reindex(sites)
    C = site_clim[BIOS].to_numpy(float)            # (nS, 19) garden climate
    c1 = C[:, 0]

    # standardize climates (founder pool for origin; site pool for garden) for match terms
    zo1 = (o1 - np.nanmean(o1)) / np.nanstd(o1)
    zc1 = (c1 - np.nanmean(c1)) / np.nanstd(c1)
    # multivariate: standardize each bioclim by founder-origin sd, project sites into same scale
    Omu, Osd = np.nanmean(O, 0), np.nanstd(O, 0) + 1e-9
    Oz = (O - Omu) / Osd; Cz = (C - Omu) / Osd
    # bio1 match and multivariate match
    M1 = -(zo1[None, :] - zc1[:, None]) ** 2                          # (nS,nF) bio1 home-field
    MM = -((Cz[:, None, :] - Oz[None, :, :]) ** 2).sum(2) / O.shape[1]  # (nS,nF) multivar

    # ---- build stacked obs (site,founder) ----
    ww = np.repeat(np.arange(nS), nF)
    ff = np.tile(np.arange(nF), nS)
    y = S.reshape(-1).astype(float)
    ok = np.isfinite(y) & np.isfinite(M1.reshape(-1))
    ww, ff, y = ww[ok], ff[ok], y[ok]
    m1 = M1.reshape(-1)[ok]; mm = MM.reshape(-1)[ok]
    y = (y - y.mean()) / y.std()                  # standardize trait
    m1 = (m1 - m1.mean()) / m1.std(); mm = (mm - mm.mean()) / mm.std()
    n = len(y)

    # eigendecomp of K once; expand to obs space: Z K Z' = (Zf K Zf')
    lamK, Vk = np.linalg.eigh(K)
    lamK = np.clip(lamK, 0, None)
    # obs-space loading: U_obs[i,:] = Vk[ff[i],:]; eigenvalues = lamK ; but each founder appears
    # in multiple sites -> U not orthonormal. Build U = (rows = obs) Vk[ff], d = lamK.
    Uobs = Vk[ff, :]                              # (n, nF)
    d = lamK.copy()

    # site fixed-effect design
    site_dum = np.zeros((n, nS)); site_dum[np.arange(n), ww] = 1.0

    def fit(match_vec):
        cols = [site_dum] + ([match_vec[:, None]] if match_vec is not None else [])
        X = np.hstack(cols)
        beta, sg2, se2, cov, Vir = reml_fit(y, X, Uobs, d)
        return beta, sg2, se2, cov, X, Vir

    # ---- in-sample variance partition + gamma test (bio1 and multivariate) ----
    res = {}
    for tag, mv in [("bio1", m1), ("multivar", mm)]:
        b0, sg0, se0, _, _, _ = fit(None)
        b1, sg1, se1, cov1, X1, _ = fit(mv)
        gamma = b1[-1]; gse = np.sqrt(cov1[-1, -1]); zg = gamma / gse
        h2 = sg1 / (sg1 + se1)
        res[tag] = dict(gamma=float(gamma), gamma_se=float(gse), gamma_z=float(zg),
                        gamma_p=float(2 * stats.norm.sf(abs(zg))),
                        h2_kinship=float(h2), sigma_g2=float(sg1), sigma_e2=float(se1),
                        h2_kinship_nomatch=float(sg0 / (sg0 + se0)))
        print(f"[{tag}] kinship h2={h2:.2f} | gamma(match)={gamma:.3f} (z={zg:.2f}, p={res[tag]['gamma_p']:.1e})")

    # ---- LEAVE-ONE-SITE-OUT CV : kinship-only vs kinship+match ----
    # Model A (kinship only): BLUP generalist value g_f -> ONE ranking for every held-out garden.
    # Model B (+match): g_f + gamma*match(o_f, c_heldout) -> re-ranks founders per garden.
    # founder BLUP is the EXACT g = K Z_f' V_unit^-1 r (sg^2 cancels), aggregated obs->founder.
    def loso(match_full):                          # match_full: (nS,nF)
        rA, rB = [], []
        for si in range(nS):
            tr = ww != si; te = ww == si
            if te.sum() < 5 or tr.sum() < 50:
                continue
            ytr = y[tr]; Btr = Uobs[tr]; ftr = ff[tr]
            keepc = [c for c in range(nS) if c != si]       # training site dummies (drop held-out)
            Xs_tr = site_dum[tr][:, keepc]
            mtr = match_full.reshape(-1)[ok][tr]
            mmu, msd = mtr.mean(), mtr.std() + 1e-9; mtr_s = (mtr - mmu) / msd
            # Model A: site means + kinship
            bA, sgA, seA, _, VirA = reml_fit(ytr, Xs_tr, Btr, d)
            aA = np.zeros(nF); np.add.at(aA, ftr, VirA); gA = K @ aA
            # Model B: + match
            XB = np.hstack([Xs_tr, mtr_s[:, None]])
            bB, sgB, seB, _, VirB = reml_fit(ytr, XB, Btr, d)
            aB = np.zeros(nF); np.add.at(aB, ftr, VirB); gB = K @ aB
            # predict held-out garden (site mean absorbed -> RANK comparison via correlation)
            yte = y[te]; fte = ff[te]
            predA = gA[fte]
            mte = (match_full[si][fte] - mmu) / msd          # standardize on TRAIN stats (no leak)
            predB = gB[fte] + bB[-1] * mte
            if np.std(yte) > 0:
                rA.append(stats.pearsonr(predA, yte)[0])
                rB.append(stats.pearsonr(predB, yte)[0])
        return np.array(rA), np.array(rB)

    cv = {}
    cv_folds = {}
    for tag, MF in [("bio1", M1), ("multivar", MM)]:
        rA, rB = loso(MF)
        cv_folds[f"{tag}_rA"] = rA; cv_folds[f"{tag}_rB"] = rB
        dr = rB - rA
        tval, pval = stats.ttest_rel(rB, rA) if len(rA) > 2 else (np.nan, np.nan)
        cv[tag] = dict(n_sites=int(len(rA)), r_kinship=float(np.mean(rA)),
                       r_kinship_match=float(np.mean(rB)), delta=float(np.mean(dr)),
                       paired_t=float(tval), paired_p=float(pval))
        print(f"[CV {tag}] kinship r={np.mean(rA):.3f} -> +match r={np.mean(rB):.3f} "
              f"(Delta={np.mean(dr):+.3f}, paired p={pval:.3f})")

    # ---- permutation null on the climate gain (shuffle founder origin labels) ----
    # Tests whether REAL origin climate beats RANDOM origin labels at re-ranking a held-out garden,
    # holding the whole CV structure fixed. Better powered than the across-fold paired t (n=17 folds
    # only). Done for both bio1 and the multivariate climate.
    rng = np.random.RandomState(0)
    NPERM = int(os.environ.get("NPERM", 200))
    perm_out, nulls = {}, {}
    for tag, MF in [("bio1", M1), ("multivar", MM)]:
        null = []
        for _ in range(NPERM):
            perm = rng.permutation(nF)
            rA, rB = loso(MF[:, perm])
            if len(rA):
                null.append(np.mean(rB - rA))
        null = np.array(null); nulls[tag] = null
        obs_delta = cv[tag]["delta"]
        perm_p = float((np.sum(null >= obs_delta) + 1) / (len(null) + 1))
        perm_out[tag] = dict(obs_delta=obs_delta, null_mean=float(null.mean()),
                             null_sd=float(null.std()), perm_p=perm_p, n_perm=int(len(null)))
        print(f"[perm {tag}] climate-gain Delta={obs_delta:+.3f}, null mean={null.mean():+.3f}, perm p={perm_p:.3f}")

    out = dict(n_sites=nS, n_founders=nF, sites=[int(s) for s in sites],
               variance_partition=res, cross_validation=cv,
               permutation=perm_out, permutation_bio1=perm_out["bio1"])  # _bio1 alias = back-compat
    json.dump(out, open(f"{H}/predict_ecotype_performance.json", "w"), indent=2)
    # per-fold CV arrays + permutation nulls -> npy so the notebook fully re-renders in `basic` env
    np.savez(f"{H}/predict_ecotype_performance_arrays.npz",
             null=nulls["bio1"], null_multivar=nulls["multivar"], **cv_folds)

    # founder-level table: generalist value (mean perf) + origin bio1 + mean match
    gen_val = np.nanmean(np.where(np.isfinite(S), S, np.nan), 0)
    tab = pd.DataFrame({"founder": founders, "origin_bio1": o1,
                        "generalist_value": gen_val})
    tab.to_csv(f"{H}/predict_ecotype_performance.csv", index=False)

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    # (1) generalist value vs origin bio1
    ax[0].scatter(o1, gen_val, s=14, c="#34495e", alpha=.6, edgecolors="none")
    rr = stats.pearsonr(o1[np.isfinite(gen_val)], gen_val[np.isfinite(gen_val)])[0]
    ax[0].set_xlabel("founder origin bio1 (home °C)"); ax[0].set_ylabel("generalist value (mean selection)")
    ax[0].set_title(f"who wins overall vs home climate  r={rr:.2f}", fontsize=10, loc="left")
    # (2) CV bars
    tags = ["bio1", "multivar"]; xb = np.arange(len(tags))
    ax[1].bar(xb - 0.2, [cv[t]["r_kinship"] for t in tags], 0.4, label="kinship only", color="#95a5a6")
    ax[1].bar(xb + 0.2, [cv[t]["r_kinship_match"] for t in tags], 0.4, label="kinship + climate-match", color="#c0392b")
    ax[1].set_xticks(xb); ax[1].set_xticklabels(tags)
    ax[1].set_ylabel("LOSO CV predictive r (held-out garden)")
    ax[1].set_title("does climate add power for a NEW garden?", fontsize=10, loc="left")
    ax[1].legend(fontsize=8, frameon=False); ax[1].axhline(0, color="k", lw=.5)
    # (3) permutation null (multivariate climate = the stronger signal)
    nb = nulls["multivar"]; ob = perm_out["multivar"]["obs_delta"]; pp = perm_out["multivar"]["perm_p"]
    ax[2].hist(nb, bins=30, color="#bdc3c7", edgecolor="none")
    ax[2].axvline(ob, color="#c0392b", lw=2, label=f"observed Δ={ob:+.3f}\nperm p={pp:.3f}")
    ax[2].set_xlabel("climate gain Δr (kinship+match − kinship)"); ax[2].set_ylabel("permutations")
    ax[2].set_title("origin-label permutation null (multivar)", fontsize=10, loc="left")
    ax[2].legend(fontsize=8, frameon=False)
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Predicting ecotype performance: 'who you are' (kinship) vs climate-match  "
                 f"({nS} gardens, {nF} founders)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(f"{H}/predict_ecotype_performance.png", dpi=150)
    print(f"\n[done] {H}/predict_ecotype_performance.{{csv,json,png}}")


if __name__ == "__main__":
    main()
