#!/usr/bin/env python
"""Pipeline B (PIPELINE_B_POOLED_MODEL.md) on the haplotype-cluster trajectories.

Per trajectory (haplotype unit), across sites:
  1. pool plots -> flower-weighted site frequency p_bar[t,g]  (t = gen 0/1/2/3)
  2. point variance V[t,g] = max(among-plot drift var of logit, binomial sampling floor)/n_plots
     gen0 (SEEDMIX founding) = among-rep var of logit(p0)/n_reps
  3. within-site WLS slope s_g of logit(p_bar) vs t, weights 1/V  -> s_g, SE(s_g)
  4. between-site WLS s_g ~ bio1 (standardized), IV weight 1/(SE^2 + tau^2) [DL];
     site-permutation null on bio1 -> two-sided + one-sided(up-in-warm) p per trajectory.

Vectorized over all trajectories (last axis). Outputs hap_gea_results.csv + GIF/Manhattan PNG.
Run in kmate env, >=16G.
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

GW = "results/grenenet_gea/gen9_window"
POOL = f"{GW}/hap/pools"
EPS = 1e-3
NPERM = 2000
RNG = np.random.default_rng(0)

def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))

# ---- load p0 (gen0 SEEDMIX hap freqs: mean + among-rep var of logit) ----
p0 = np.load(f"{GW}/hap/p0_seedmix.npy")            # [n_traj] mean founding hap freq
v0 = np.load(f"{GW}/hap/v0_seedmix.npy")            # [n_traj] among-rep var of logit(p0)/n_reps
N = p0.shape[0]

# ---- load evolved-gen pool matrices + meta ----
gens = [1, 2, 3]
pools = {g: np.load(f"{POOL}/pool_gen{g}_hap.npy") for g in gens}      # [n_pools_g x N]
meta = {g: pd.read_csv(f"{POOL}/pool_gen{g}.meta.csv") for g in gens}
sites = sorted(set().union(*[set(meta[g].site) for g in gens]))
# per-site bio1 (one value per site)
site_bio1 = {}
for g in gens:
    for _, r in meta[g].iterrows():
        site_bio1[int(r.site)] = float(r.bio1)

# ---- Steps 1-2: site freq p_bar[t] and variance V[t] per site ----
# accumulate, per site, the trajectory arrays across t=0..3
T = [0, 1, 2, 3]
site_data = {}                       # site -> (t_list, Pbar[k x N], V[k x N])
sig_pool = None                      # pooled drift var for shrinkage (computed below)
# first pass: gather per-site,per-gen pooled freq + drift var
tmp = {s: {} for s in sites}
drift_all = []
for g in gens:
    m = meta[g]; P = pools[g]
    for s in sites:
        idx = np.where(m.site.values == s)[0]
        if len(idx) == 0:
            continue
        fl = m.total_flowers.values[idx][:, None]            # [n_plots x 1]
        cov = m.mean_coverage.values[idx][:, None]
        pp = P[idx]                                          # [n_plots x N]
        pbar = (fl * pp).sum(0) / fl.sum()                   # [N] flower-weighted site freq
        lp = logit(pp)                                       # [n_plots x N]
        if len(idx) >= 2:
            drift = lp.var(0, ddof=1)                        # among-plot drift var of logit
            drift_all.append(drift)
        else:
            drift = np.full(N, np.nan)
        # binomial sampling floor (per-plot, averaged): Neff = 1/(1/(2 fl)+1/cov)
        neff = 1.0 / (1.0 / (2 * fl) + 1.0 / np.maximum(cov, 1e-6))   # [n_plots x 1]
        vsamp = (1.0 / (neff * np.clip(pbar, EPS, 1 - EPS) * (1 - np.clip(pbar, EPS, 1 - EPS)))).mean(0)
        tmp[s][g] = dict(pbar=pbar, drift=drift, vsamp=vsamp, n=len(idx))
sig_pool = np.nanmedian(np.vstack(drift_all), axis=0) if drift_all else np.full(N, 0.1)

# assemble per-site (t, Pbar, V), gen0 anchor = p0 for every site
for s in sites:
    ts, Pl, Vl = [0], [logit(p0)], [v0.copy()]
    for g in gens:
        if g not in tmp[s]:
            continue
        d = tmp[s][g]; n = d["n"]
        drift = d["drift"]
        drift_sh = np.where(np.isnan(drift), sig_pool, ((n - 1) * np.nan_to_num(drift) + 4 * sig_pool) / ((n - 1) + 4))
        V = np.maximum(drift_sh, d["vsamp"]) / max(n, 1)
        ts.append(g); Pl.append(logit(d["pbar"])); Vl.append(V)
    site_data[s] = (np.array(ts), np.vstack(Pl), np.vstack(Vl))   # [k], [k x N], [k x N]

# ---- Step 3: per-site weighted slope s_g + SE (vectorized over N) ----
def wls_slope(tarr, Y, V):
    """tarr [k]; Y,V [k x N] -> slope[N], se[N] of weighted lstsq Y ~ t."""
    w = 1.0 / V                                              # [k x N]
    tw = tarr[:, None]
    sw = w.sum(0)
    tbar = (w * tw).sum(0) / sw
    Sxx = (w * (tw - tbar) ** 2).sum(0)
    slope = (w * (tw - tbar) * Y).sum(0) / Sxx
    se = np.sqrt(1.0 / Sxx)
    return slope, se

S = np.full((len(sites), N), np.nan); SE = np.full((len(sites), N), np.nan)
for i, s in enumerate(sites):
    ts, Y, V = site_data[s]
    if len(ts) >= 2:
        S[i], SE[i] = wls_slope(ts, Y, V)
clim = np.array([site_bio1[s] for s in sites], float)
clim = (clim - clim.mean()) / clim.std()

# ---- Step 4: s_g ~ bio1, IV-weighted with DL tau^2; permutation null ----
def meta_beta(S, SE, x):
    """S,SE [G x N]; x [G] -> beta1[N] (IV-weighted, DL tau^2)."""
    G = S.shape[0]
    w0 = 1.0 / SE ** 2
    # DL tau^2 from a fixed-effect intercept+slope fit residuals
    xc = x[:, None]
    sw = w0.sum(0); xbar = (w0 * xc).sum(0) / sw
    Sxx = (w0 * (xc - xbar) ** 2).sum(0)
    b1 = (w0 * (xc - xbar) * S).sum(0) / Sxx
    b0 = (w0 * S).sum(0) / sw - b1 * xbar
    resid = S - (b0 + b1 * xc)
    Q = (w0 * resid ** 2).sum(0)
    tau2 = np.maximum(0.0, (Q - (G - 2)) / (sw - (w0 ** 2).sum(0) / sw))
    w = 1.0 / (SE ** 2 + tau2)
    swt = w.sum(0); xbart = (w * xc).sum(0) / swt
    Sxxt = (w * (xc - xbart) ** 2).sum(0)
    return (w * (xc - xbart) * S).sum(0) / Sxxt

valid = np.isfinite(S).sum(0) >= 5            # need >=5 sites with a slope
beta = meta_beta(np.where(np.isfinite(S), S, 0), np.where(np.isfinite(SE), SE, 1e6), clim)
cnt_two = np.zeros(N); cnt_up = np.zeros(N)
Sf = np.where(np.isfinite(S), S, 0); SEf = np.where(np.isfinite(SE), SE, 1e6)
for _ in range(NPERM):
    xp = clim[RNG.permutation(len(sites))]
    bp = meta_beta(Sf, SEf, xp)
    cnt_two += (np.abs(bp) >= np.abs(beta))
    cnt_up += (bp >= beta)
p_two = (cnt_two + 1) / (NPERM + 1)
p_up = (cnt_up + 1) / (NPERM + 1)

reg = pd.read_csv(f"{GW}/hap/registry.csv")
reg["beta1_climate"] = beta
reg["p_perm_two"] = p_two
reg["p_perm_up"] = p_up
reg["n_sites"] = np.isfinite(S).sum(0)
reg["testable"] = valid
out = reg[reg.testable].copy()
out.to_csv(f"{GW}/hap/hap_gea_results.csv", index=False)

# GIF (genomic inflation) from two-sided perm p
from scipy import stats
chi = stats.chi2.isf(out.p_perm_two.clip(1e-6, 1), 1)
gif = np.median(chi) / stats.chi2.ppf(0.5, 1)
print(f"testable trajectories: {out.shape[0]:,} / {N:,}")
print(f"GIF (lambda) = {gif:.3f}")
print(f"raw p<0.05: {100*(out.p_perm_two<0.05).mean():.1f}%  | one-sided up-in-warm p<0.05: {100*(out.p_perm_up<0.05).mean():.1f}%")
print(f"top hits (smallest two-sided perm p):")
print(out.sort_values('p_perm_two')[['chrom','start','end','founding_freq','beta1_climate','p_perm_two','p_perm_up','n_sites']].head(12).to_string(index=False))
print(f"\nwrote {GW}/hap/hap_gea_results.csv")
