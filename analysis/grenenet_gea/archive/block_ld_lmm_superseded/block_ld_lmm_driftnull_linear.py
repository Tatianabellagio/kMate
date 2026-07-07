#!/usr/bin/env python
"""WF-drift-through-projection null for the LINEAR + sampling-floor block LD-LMM (default site 4).

Same selfing-aware generative null as block_ld_lmm_driftnull.py, but the per-sim pipeline is the
LINEAR-scale + binomial-floor model (block_ld_lmm_linear_sampvar.py) and the test is ONE-SIDED for
the target "rose MORE than founder-linkage predicts" (resid>0). Reports:
  * upper-tail maxT FWER threshold (95th pct of max_b z_b over sims) for rose-MORE, over the
    BOUNDARY-GUARDED block set (panel_freq in [GLO,GHI]); per-block p_pos_fwer;
  * two-sided max|z| threshold too (for reference).
REML-refit per sim (self-consistent scaling). Writes site<ID>_block_ld_lmm_driftnull_linear.csv
+ _meta.json. Env: kmate. SITE, NSIM via env.
"""
import os, sys, json, glob
import numpy as np
import pandas as pd
from scipy import stats, optimize
from scipy.linalg import cho_factor, cho_solve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from founder_genotype import build_genotype

H_ = "results/grenenet_gea/hapfreq"
WIN = "results/grenenet_kmate_window"
SEED = "results/grenenet_kmate_window_seedmix"
SITE = int(os.environ.get("SITE", 4))
NSIM = int(os.environ.get("NSIM", 2000))
GLO, GHI = 0.10, 0.90                                   # boundary guard for the target set
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
Tg = np.array([0.0, 1.0, 2.0, 3.0]); coef = (Tg - 1.5)
rng = np.random.RandomState(0)


def bvar(f, N):
    f = np.clip(f, 1e-6, 1 - 1e-6)
    return f * (1 - f) / np.clip(N, 1.0, None)


def genome_h(samp, base):
    gs = []
    for ch in CHROMS:
        f = f"{base}/{samp}_{ch}.h_blocks_per_chrom.npz"
        if not os.path.exists(f):
            return None
        gs.append(np.load(f, allow_pickle=True)[f"{ch}_global_h"].astype(np.float64))
    return np.mean(gs, 0)


def main():
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
    pfreq = fr.loc[ki, "panel_freq"].to_numpy()
    guard = (pfreq >= GLO) & (pfreq <= GHI)                # target block set
    Graw = G[:, ki].astype(np.float64); nF = Graw.shape[0]

    win_founders = np.load(glob.glob(f"{SEED}/*_Chr1.h_blocks_per_chrom.npz")[0],
                           allow_pickle=True)["founders"].astype(str)
    assert np.array_equal(np.asarray(founders).astype(str), win_founders), "founder order mismatch!"

    # founder trajectories + Ne
    seeds = sorted({p.split("/")[-1].split("_Chr")[0]
                    for p in glob.glob(f"{SEED}/*_Chr1.h_blocks_per_chrom.npz")})
    P0 = np.vstack([genome_h(s, SEED) for s in seeds])
    p0f = P0.mean(0); v0f = P0.var(0, ddof=1)
    pt = lib.pool_table()
    s = pt[pt.site == SITE].copy()
    s = s[s.sampleid.astype(str).apply(lambda x: os.path.exists(f"{WIN}/{x}_Chr1.h_blocks_per_chrom.npz"))]
    cellH, cellN = {}, {}
    for (gen, plot), g in s.groupby(["generation", "plot"]):
        hs, ws = [], []
        for _, r in g.iterrows():
            h = genome_h(str(r.sampleid), WIN)
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
    flowers = np.ones((4, n))
    for ti, gen in enumerate([1, 2, 3], start=1):
        Hbar[ti] = np.mean([cellH[(gen, pl)] for pl in plots], 0)
        for j, pl in enumerate(plots):
            flowers[ti, j] = cellN[(gen, pl)]
    Nmed = np.median(flowers[1:].ravel()); flowers[0] = Nmed
    H3 = np.vstack([cellH[(3, pl)] for pl in plots]); freq3 = H3.mean(0)
    res = (freq3 > 0.02) & (freq3 < 0.4)
    A_ = float(np.clip(np.median(H3[:, res].var(0, ddof=1) / (freq3[res] * (1 - freq3[res]))), 1e-4, 0.99))
    Ne = float(np.clip(1.0 / (1.0 - (1.0 - A_) ** (1.0 / 3.0)), 5, 1e5))

    A = Graw.copy(); pf = A.mean(0); A = (A - pf) / np.sqrt(pf * (1 - pf) + 1e-9)
    A32 = A.astype(np.float32); one = np.ones(M); _warm = [None]

    def fit_z(freq, refit=True):
        """freq: (n,4,M) RAW block hap freqs -> block residual z (LINEAR scale, binomial floor)."""
        slopes = np.zeros((n, M)); evo = np.zeros(M)
        for j in range(n):
            y = np.clip(freq[j], 0.0, 1.0)
            slopes[j] = (coef[:, None] * y).sum(0) / 5.0
            for ti in (1, 2, 3):
                evo += (coef[ti] / 5.0) ** 2 * bvar(freq[j, ti], 2 * flowers[ti, j])
        evo /= n ** 2
        v0_term = (coef[0] / 5.0) ** 2 * bvar(freq[:, 0, :].mean(0), 2 * Nmed)
        sb = slopes.mean(0); sb = sb - sb.mean()
        rep = slopes.std(0, ddof=1) ** 2 / n
        D0 = v0_term + np.maximum(rep, evo); D0 = np.clip(D0, np.median(D0) * 1e-3, None)

        def pieces(tau, se2):
            Winv = 1.0 / (se2 + D0); Wv32 = Winv.astype(np.float32)
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
            _warm[0] = np.array([tau, se2])
        Vinv, _ = pieces(tau, se2)
        mu = (one @ Vinv(sb)) / (one @ Vinv(one)); r = sb - mu
        g_b = (tau / nF) * (A.T @ (A @ Vinv(r)))
        return (r - g_b) / np.sqrt(se2 + D0), float(tau), float(se2)

    # observed
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
    z_obs, tau0, se20 = fit_z(obs)
    print(f"  observed REML refit: tau={tau0:.5f}, se2={se20:.5f}; obs max z (guarded)={z_obs[guard].max():.2f}")

    NULLMODE = os.environ.get("NULLMODE", "walk")          # "walk" (drift trajectory) | "indep" (old)

    def sim_freq():
        Hs = np.empty((n, 4, nF))
        # gen0: SHARED across plots (mirrors the single seedmix p0; per-sim founding shift only)
        Hs[:, 0, :] = np.clip(Hbar[0] + rng.normal(0, np.sqrt(v0f)), 0, 1)[None, :]
        if NULLMODE == "walk":
            # WF drift TRAJECTORY: f_t = WFstep(f_{t-1}) + deterministic winning increment (tracks
            # Hbar in mean, accumulates drift -> autocorrelated, Ne-consistent with the gen3 estimator)
            for ti in (1, 2, 3):
                prev = np.clip(Hs[:, ti - 1, :], 0, 1)
                drifted = rng.binomial(int(Ne), prev) / Ne
                Hs[:, ti, :] = np.clip(drifted + (Hbar[ti] - Hbar[ti - 1])[None, :], 0, 1)
        else:
            for ti in (1, 2, 3):
                Hs[:, ti, :] = rng.binomial(int(Ne), np.clip(Hbar[ti], 0, 1), size=(n, nF)) / Ne
        Hs /= np.clip(Hs.sum(-1, keepdims=True), 1e-12, None)
        return np.einsum("ptf,fb->ptb", Hs, Graw)

    maxz = np.empty(NSIM); maxabs = np.empty(NSIM)
    ge_pos = np.zeros(M, np.int64); ge_abs = np.zeros(M, np.int64); pool = []
    for k in range(NSIM):
        zk, _, _ = fit_z(sim_freq())
        zg = zk[guard]
        mxz = float(zg.max()); mxa = float(np.abs(zg).max())
        maxz[k] = mxz; maxabs[k] = mxa
        ge_pos += (mxz >= z_obs); ge_abs += (mxa >= np.abs(z_obs))
        pool.append(zg[rng.randint(0, guard.sum(), size=128)])
        if (k + 1) % 200 == 0:
            print(f"  ...{k+1}/{NSIM} sims  (rose-more maxT 95% = {np.quantile(maxz[:k+1],0.95):.2f})")
    pool = np.concatenate(pool)

    thr_pos = float(np.quantile(maxz, 0.95)); thr_abs = float(np.quantile(maxabs, 0.95))
    p_pos_fwer = (1 + ge_pos) / (1 + NSIM)
    ps = np.sort(pool)
    p_pos_marg = np.clip(1.0 - np.searchsorted(ps, z_obs, side="right") / len(ps), 1.0 / len(ps), 1.0)

    d = fr.loc[ki, ["chrom", "unit_start", "unit_end", "panel_freq"]].copy().reset_index(drop=True)
    d["gid"] = ki                                          # unique global haplotype id (join key; units hold k-1 haps)
    d["z"] = z_obs; d["guard"] = guard
    d["p_pos_marg"] = p_pos_marg; d["p_pos_fwer"] = p_pos_fwer
    # one-sided BH over guarded blocks only (positional fill -> no unit-key merge duplication)
    qfull = np.full(M, np.nan)
    gi = np.where(guard)[0]; ppm = p_pos_marg[gi]; mg = len(gi); o = ppm.argsort()
    qq = np.empty(mg); qq[o] = np.minimum.accumulate((ppm[o] * mg / (np.arange(mg) + 1))[::-1])[::-1]
    qfull[gi] = np.clip(qq, 0, 1)
    d["q_pos_drift"] = qfull
    d.sort_values("p_pos_marg").to_csv(f"{H_}/site{SITE}_block_ld_lmm_driftnull_linear.csv", index=False)

    nfw = int((d.guard & (d.p_pos_fwer < 0.05)).sum())
    meta = dict(site=SITE, n_plot=n, M=int(M), n_guard=int(guard.sum()), NSIM=NSIM, Ne=float(Ne),
                nullmode=NULLMODE,
                tau=tau0, se2=se20, guard=[GLO, GHI],
                rose_more_maxT_fwer05=thr_pos, twosided_maxT_fwer05=thr_abs,
                obs_max_z_guarded=float(z_obs[guard].max()),
                pool_z_q={q_: float(np.quantile(pool, q_)) for q_ in [0.5, 0.95, 0.99, 0.999]},
                n_rose_more_fwer05=nfw,
                n_rose_more_qdrift05=int((d.q_pos_drift < 0.05).sum()))
    json.dump(meta, open(f"{H_}/site{SITE}_block_ld_lmm_driftnull_linear_meta.json", "w"), indent=2)

    print(f"\nsite {SITE} LINEAR drift null: {M:,} blocks ({int(guard.sum())} guarded [{GLO},{GHI}]), "
          f"Ne~{Ne:.0f}, {NSIM} sims")
    print(f"  rose-MORE one-sided maxT FWER 5% threshold:  z > {thr_pos:.2f}   (obs max guarded z = {z_obs[guard].max():.2f})")
    print(f"  rose-more hits: FWER p<0.05 = {nfw} | drift-BH q<0.05 = {int((d.q_pos_drift<0.05).sum())}")
    print("\n  TOP rose-MORE blocks (guarded, by drift one-sided p):")
    gg = d[d.guard].sort_values("p_pos_marg").head(14)
    print(gg[["chrom", "unit_start", "unit_end", "panel_freq", "z", "p_pos_marg", "p_pos_fwer", "q_pos_drift"]].to_string(index=False))
    print(f"\n[done] {H_}/site{SITE}_block_ld_lmm_driftnull_linear.csv + _meta.json")


if __name__ == "__main__":
    main()
