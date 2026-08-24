#!/usr/bin/env python
"""WF-drift-through-projection NULL + FWER threshold for the sampling-variance block LD-LMM.

The sign-flip permutation is invalid here (it annihilates the founder-level selection that must
be PRESERVED under H0 -> anti-conservative). The correct null is GENERATIVE and exploits that,
under near-complete selfing, whole ecotype genomes segregate as units: H0 ("no BLOCK-specific
selection") == "every block is a pure projection of its founders". So we simulate the founder
frequencies and project -- which captures the cross-block residual correlation the parametric
null misses, and gets SELFING for free (no selfing parameter):

  * SELFING's drift inflation is absorbed into the pooled Ne fit EMPIRICALLY from the observed
    among-plot variance of gen-3 frequencies (same estimator as ecotype_selection_site.py).
  * SELFING's whole-genome segregation IS the null: we draw founder-level (231-dim) drift +
    flower-census sampling noise around the REAL winning trajectory H_bar (preserves the winners,
    hence the boundary curvature & frequency levels), then PROJECT f_block = sum_{f in b} h_f.
    Sampling a whole genome moves all its blocks together -> correct correlated cross-block noise.
    No block-specific term is injected -> any residual is pure null (drift+sampling+curvature+proj).

Per simulated "site": for each plot j and gen t, H_sim[j,t,:] = binom(Ne, H_bar[t,:])/Ne
(t=0 also gets seedmix founding error), renormalised to the simplex; project via the founder x
haplotype membership G; run the IDENTICAL sampling-variance Stage-1/Stage-2 at the FITTED
(tau, sigma_e^2); record block residual z. maxT over NSIM sims -> FWER threshold; per-block
FWER p and a pooled marginal null. Writes site<ID>_block_ld_lmm_driftnull.csv + _meta.json.
Env: kmate. SITE, NSIM via env.
"""
import os, sys, json, glob
import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import cho_factor, cho_solve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from founder_genotype import build_genotype

H_ = "results/grenenet_gea/hapfreq"
WIN = "results/grenenet_kmate_window"
SEED = "results/grenenet_kmate_window_seedmix"
SITE = int(os.environ.get("SITE", 4))
NSIM = int(os.environ.get("NSIM", 500))
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
EPS = 1e-3
Tg = np.array([0.0, 1.0, 2.0, 3.0])
coef = (Tg - 1.5)                                       # sum sq = 5
rng = np.random.RandomState(0)


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def logit_var(f, Ngam):
    f = np.clip(f, EPS, 1 - EPS)
    return 1.0 / (np.clip(Ngam, 1.0, None) * f * (1 - f))


def genome_h(samp, base, founders_ref):
    gs = []
    for ch in CHROMS:
        f = f"{base}/{samp}_{ch}.h_blocks_per_chrom.npz"
        if not os.path.exists(f):
            return None
        gs.append(np.load(f, allow_pickle=True)[f"{ch}_global_h"].astype(np.float64))
    return np.mean(gs, 0)


def main():
    # ---- membership / projection (founder x testable-haplotype), aligned founder order ----
    G, founders, reg = build_genotype()
    mat = np.load(f"{H_}/hapfreq_matrix.npy")
    samples = [l.strip() for l in open(f"{H_}/hapfreq_samples.txt") if l.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    p0_hap = np.load(f"{H_}/hapfreq_p0_seedmix.npy")
    fr = pd.read_csv(f"{H_}/hapfreq_registry.csv")
    nhap = mat.shape[1]
    fr["unit"] = fr.chrom + ":" + fr.unit_start.astype(str) + "-" + fr.unit_end.astype(str)
    testable = (fr.covered & fr.panel_freq.between(0.05, 0.95)).to_numpy()
    keep = np.zeros(nhap, bool)
    for u, grp in fr[testable].groupby("unit"):
        kept = grp.index.drop(grp.panel_freq.idxmax()) if len(grp) > 1 else grp.index
        keep[kept.to_numpy()] = True
    ki = np.where(keep)[0]; M = len(ki)
    Graw = G[:, ki].astype(np.float64)                  # (nF, M) 0/1 membership -> projection
    nF = Graw.shape[0]

    # founder order check (membership vs window h)
    win_founders = np.load(glob.glob(f"{SEED}/*_Chr1.h_blocks_per_chrom.npz")[0],
                           allow_pickle=True)["founders"].astype(str)
    assert np.array_equal(np.asarray(founders).astype(str), win_founders), "founder order mismatch!"

    # ---- founder trajectories H[j,t,f] + Ne (same as ecotype_selection_site.py) ----
    seeds = sorted({p.split("/")[-1].split("_Chr")[0]
                    for p in glob.glob(f"{SEED}/*_Chr1.h_blocks_per_chrom.npz")})
    P0 = np.vstack([genome_h(s, SEED, founders) for s in seeds])
    p0f = P0.mean(0); v0f = P0.var(0, ddof=1)           # founder founding mean + error
    pt = lib.pool_table()
    s = pt[pt.site == SITE].copy()
    s = s[s.sampleid.astype(str).apply(
        lambda x: os.path.exists(f"{WIN}/{x}_Chr1.h_blocks_per_chrom.npz"))]
    cellH, cellN = {}, {}
    for (gen, plot), g in s.groupby(["generation", "plot"]):
        hs, ws = [], []
        for _, r in g.iterrows():
            h = genome_h(str(r.sampleid), WIN, founders)
            if h is not None:
                w = r.flowerscollected if np.isfinite(r.flowerscollected) and r.flowerscollected > 0 else 1.0
                hs.append(h); ws.append(w)
        if hs:
            ws = np.array(ws)
            cellH[(int(gen), int(plot))] = (np.vstack(hs) * ws[:, None]).sum(0) / ws.sum()
            cellN[(int(gen), int(plot))] = float(ws.sum())
    present = {}
    for (gen, plot) in cellH:
        present.setdefault(plot, set()).add(gen)
    plots = sorted(pl for pl, gs in present.items() if {1, 2, 3} <= gs)
    n = len(plots)
    Hbar = np.zeros((4, nF)); Hbar[0] = p0f
    flowers = np.ones((4, n))                            # founding flowers ~ Nmed proxy filled below
    for ti, gen in enumerate([1, 2, 3], start=1):
        Hbar[ti] = np.mean([cellH[(gen, pl)] for pl in plots], 0)
        for j, pl in enumerate(plots):
            flowers[ti, j] = cellN[(gen, pl)]
    Nmed = np.median(flowers[1:].ravel())
    flowers[0] = Nmed
    # pooled Ne from among-plot variance of well-resolved founders at gen 3
    H3 = np.vstack([cellH[(3, pl)] for pl in plots])
    freq3 = H3.mean(0)
    res = (freq3 > 0.02) & (freq3 < 0.4)
    A_ = float(np.clip(np.median(H3[:, res].var(0, ddof=1) / (freq3[res] * (1 - freq3[res]))), 1e-4, 0.99))
    Ne = float(np.clip(1.0 / (1.0 - (1.0 - A_) ** (1.0 / 3.0)), 5, 1e5))

    # ---- LD-LMM machinery: REML-REFIT (tau, se2) on EACH dataset (observed AND every sim),
    #      so processing is identical on both sides and the null z is self-consistently scaled. ----
    from scipy import optimize
    A = Graw.copy(); pf = A.mean(0); A = (A - pf) / np.sqrt(pf * (1 - pf) + 1e-9)
    A32 = A.astype(np.float32)                            # float32 for the hot AWA matmul (2x)
    one = np.ones(M)
    _warm = [None]                                        # warm-start (tau,se2) across sims

    def fit_z(freq, refit=True):
        """freq: (n,4,M) per-plot per-gen block hap freqs -> block residual z (REML-refit)."""
        slopes = np.zeros((n, M)); evo = np.zeros(M)
        for j in range(n):
            y = logit(freq[j])                           # (4,M)
            slopes[j] = (coef[:, None] * y).sum(0) / 5.0
            for ti in (1, 2, 3):
                evo += (coef[ti] / 5.0) ** 2 * logit_var(freq[j, ti], 2 * flowers[ti, j])
        evo /= n ** 2
        v0_term = (coef[0] / 5.0) ** 2 * logit_var(freq[:, 0, :].mean(0), 2 * Nmed)
        sb = slopes.mean(0); sb = sb - sb.mean()
        rep = slopes.std(0, ddof=1) ** 2 / n
        D0 = v0_term + np.maximum(rep, evo)
        D0 = np.clip(D0, np.median(D0) * 1e-3, None)

        def pieces(tau, se2):
            Winv = 1.0 / (se2 + D0)
            Wv32 = Winv.astype(np.float32)
            AWA = ((A32 * Wv32[None, :]) @ A32.T).astype(np.float64) / nF
            cf = cho_factor(np.eye(nF) / tau + AWA, lower=True)

            def Vinv(x):
                Wx = Winv * x
                u = cho_solve(cf, (A @ Wx) / np.sqrt(nF))
                return Wx - Winv * ((A.T @ u) / np.sqrt(nF))
            logdetV = np.sum(np.log(se2 + D0)) + np.linalg.slogdet(np.eye(nF) + tau * AWA)[1]
            return Vinv, logdetV

        def neg_reml(par):
            tau, se2 = np.exp(par)
            Vinv, logdetV = pieces(tau, se2)
            Vi1 = Vinv(one); mu = (one @ Vinv(sb)) / (one @ Vi1)
            r = sb - mu
            return 0.5 * (logdetV + r @ Vinv(r) + np.log(one @ Vi1))

        x0 = np.log(_warm[0]) if (_warm[0] is not None and refit) else np.log([np.var(sb), np.var(sb) / 2])
        res = optimize.minimize(neg_reml, x0=x0, method="Nelder-Mead",
                                options={"xatol": 1e-2, "fatol": 1e-2, "maxiter": 80})
        tau, se2 = np.exp(res.x)
        if refit and _warm[0] is None:
            _warm[0] = np.array([tau, se2])              # set warm start from the observed fit
        Vinv, _ = pieces(tau, se2)
        Vi1 = Vinv(one); mu = (one @ Vinv(sb)) / (one @ Vi1)
        r = sb - mu
        g_b = (tau / nF) * (A.T @ (A @ Vinv(r)))
        return (r - g_b) / np.sqrt(se2 + D0), float(tau), float(se2)

    resid_z = lambda freq: fit_z(freq)[0]

    # ---- OBSERVED z (recomputed identically; sanity-checked vs the sampvar csv) ----
    obs = np.zeros((n, 4, M)); obs[:, 0, :] = p0_hap[ki][None, :]
    s2 = pt[(pt.site == SITE) & pt.sampleid.astype(str).isin(sidx)].copy()
    ocell = {}
    for (gen, plot), g in s2.groupby(["generation", "plot"]):
        rows = g.sampleid.astype(str).map(sidx).to_numpy()
        w = g.flowerscollected.to_numpy(float); w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
        ocell[(int(gen), int(plot))] = ((mat[rows] * w[:, None]).sum(0) / w.sum())[ki]
    for j, pl in enumerate(plots):
        for ti, gen in enumerate([1, 2, 3], start=1):
            obs[j, ti, :] = ocell[(gen, pl)]
    z_obs, tau0, se20 = fit_z(obs); absz_obs = np.abs(z_obs)
    print(f"  observed REML refit: tau={tau0:.4f}, se2={se20:.4f}  (sampvar meta ~0.0589/0.0068); "
          f"obs max|z|={absz_obs.max():.2f}")

    # ---- DRIFT-PROJECTION NULL ----
    def sim_freq():
        Hs = np.empty((n, 4, nF))
        for ti in range(4):
            base = np.clip(Hbar[ti] + (rng.normal(0, np.sqrt(v0f)) if ti == 0 else 0.0), 0, 1)
            draw = rng.binomial(int(Ne), np.clip(base, 0, 1), size=(n, nF)) / Ne
            Hs[:, ti, :] = draw
        Hs /= np.clip(Hs.sum(-1, keepdims=True), 1e-12, None)   # back to the simplex (per plot,gen)
        return np.einsum("ptf,fb->ptb", Hs, Graw)               # project -> (n,4,M)

    maxabs = np.empty(NSIM); ge = np.zeros(M, np.int64); pool = []
    for k in range(NSIM):
        zk = np.abs(resid_z(sim_freq()))
        mk = float(zk.max()); maxabs[k] = mk
        ge += (mk >= absz_obs)
        pool.append(zk[rng.randint(0, M, size=128)])
        if (k + 1) % 100 == 0:
            print(f"  ...{k+1}/{NSIM} sims  (running maxT 95% = {np.quantile(maxabs[:k+1],0.95):.2f})")
    pool = np.concatenate(pool)

    thr05 = float(np.quantile(maxabs, 0.95)); thr10 = float(np.quantile(maxabs, 0.90))
    p_fwer = (1 + ge) / (1 + NSIM)
    ps = np.sort(pool)
    p_marg = np.clip(1.0 - np.searchsorted(ps, absz_obs, side="right") / len(ps), 1.0 / len(ps), 1.0)

    d = fr.loc[ki, ["chrom", "unit_start", "unit_end", "panel_freq"]].copy()
    d["z"] = z_obs; d["p_param"] = 2 * stats.norm.sf(absz_obs)
    d["p_marg_drift"] = p_marg; d["p_fwer"] = p_fwer
    m = len(d); o = d.p_marg_drift.values.argsort()
    q = np.empty(m); q[o] = np.minimum.accumulate((d.p_marg_drift.values[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    d["q_drift"] = np.clip(q, 0, 1)
    d.sort_values("p_marg_drift").to_csv(f"{H_}/site{SITE}_block_ld_lmm_driftnull.csv", index=False)

    meta = dict(site=SITE, n_plot=n, M=int(M), nF=int(nF), NSIM=NSIM, Ne=float(Ne),
                tau=tau0, se2=se20, maxT_thr_fwer05=thr05, maxT_thr_fwer10=thr10,
                obs_max_absz=float(absz_obs.max()),
                pool_z_q={q_: float(np.quantile(pool, q_)) for q_ in [0.5, 0.95, 0.99, 0.999]},
                n_fwer05=int((d.p_fwer < 0.05).sum()), n_qdrift05=int((d.q_drift < 0.05).sum()))
    json.dump(meta, open(f"{H_}/site{SITE}_block_ld_lmm_driftnull_meta.json", "w"), indent=2)

    print(f"\nsite {SITE}: {M:,} blocks, {n} plots, Ne~{Ne:.0f}, {NSIM} drift-projection sims")
    print(f"  observed max |z| = {absz_obs.max():.2f}")
    print(f"  drift null |z| quantiles 50/95/99/99.9 = "
          + "/".join(f"{np.quantile(pool, q_):.2f}" for q_ in [0.5, 0.95, 0.99, 0.999])
          + "   (param N(0,1): 0.67/1.96/2.58/3.29)")
    print(f"  maxT FWER threshold: 5% |z| > {thr05:.2f}   10% |z| > {thr10:.2f}")
    print(f"  hits: FWER p<0.05 = {int((d.p_fwer<0.05).sum())} | drift-BH q<0.05 = {int((d.q_drift<0.05).sum())}")
    print("\n  TOP blocks by drift-null marginal p:")
    cols = ["chrom", "unit_start", "unit_end", "panel_freq", "z", "p_param", "p_marg_drift", "p_fwer", "q_drift"]
    print(d.sort_values("p_marg_drift").head(12)[cols].to_string(index=False))
    print(f"\n[done] {H_}/site{SITE}_block_ld_lmm_driftnull.csv + _meta.json")


if __name__ == "__main__":
    main()
