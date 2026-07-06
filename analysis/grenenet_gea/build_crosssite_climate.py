#!/usr/bin/env python
"""Cross-site climate association of per-block founding->gen1 selection (the multi-site combination).

For each haploblock we have a per-site selection vector from block_ld_lmm_gen1.py (run over all 30
sites). We test, per block, whether that selection TRACKS GARDEN CLIMATE across sites, for BOTH:
  * s_b  = raw block selection (includes ecotype hitchhiking)  -> ~ecotype-level climate adaptation
  * r_b  = block-specific residual (beyond founder linkage)    -> genuinely block-specific selection

Statistic = site-weighted (w=n_plot) Pearson correlation of the block's site-vector with climate.
NULL = site-label PERMUTATION of the climate vector (shuffle the 30 garden climates): preserves the
among-block LD structure, needs no per-site drift threshold. maxT over blocks -> FWER; pooled ->
marginal p -> BH q. Run for climate PC1 (main gradient) and bio1 (annual mean temp). Compares how
much of the climate signal is ecotype-level (s) vs block-specific (r). Env: kmate.
"""
import os, sys, glob, json
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
H = "results/grenenet_gea/hapfreq"
CLIM = "/global/scratch/users/tbellg/gea_grene-net/bioclimvars_experimental_sites_era5.csv"
NPERM = int(os.environ.get("NPERM", 5000))
rng = np.random.RandomState(0)


def wcorr_setup(E_raw, w):
    """precompute weighted-centered, normalized effect rows E (blocks x sites)."""
    ebar = (E_raw * w[None, :]).sum(1, keepdims=True) / w.sum()
    Ew = np.sqrt(w)[None, :] * (E_raw - ebar)
    nrm = np.sqrt((Ew ** 2).sum(1, keepdims=True)); nrm[nrm == 0] = 1.0
    return Ew / nrm                                       # (B, S) unit weighted-centered rows


def wcorr_vec(En, c, w):
    """weighted corr of every block row (En precomputed) with climate vector c."""
    cbar = (w * c).sum() / w.sum()
    Cw = np.sqrt(w) * (c - cbar); n = np.sqrt((Cw ** 2).sum())
    return En @ (Cw / (n if n else 1.0))                 # (B,)


def run_assoc(En, clim, w, label):
    obs = wcorr_vec(En, clim, w)
    B = len(obs)
    maxabs = np.empty(NPERM); ge = np.zeros(B, np.int64); pool = []
    for k in range(NPERM):
        cp = clim[rng.permutation(len(clim))]
        r = wcorr_vec(En, cp, w)
        a = np.abs(r); mx = a.max(); maxabs[k] = mx
        ge += (mx >= np.abs(obs))
        pool.append(a[rng.randint(0, B, size=64)])
    pool = np.sort(np.concatenate(pool))
    p_fwer = (1 + ge) / (1 + NPERM)
    p_marg = np.clip(1.0 - np.searchsorted(pool, np.abs(obs), "right") / len(pool), 1.0 / len(pool), 1.0)
    o = p_marg.argsort(); B_ = B
    q = np.empty(B_); q[o] = np.minimum.accumulate((p_marg[o] * B_ / (np.arange(B_) + 1))[::-1])[::-1]
    thr = float(np.quantile(maxabs, 0.95))
    print(f"  [{label}] maxT FWER5% |r|>{thr:.3f}; FWER hits {int((p_fwer<0.05).sum())}; "
          f"BH q<0.05 {int((np.clip(q,0,1)<0.05).sum())}; max|r|={np.abs(obs).max():.3f}")
    return obs, p_marg, p_fwer, np.clip(q, 0, 1), thr


def main():
    # ---- assemble per-site effect matrices aligned by gid ----
    files = sorted(glob.glob(f"{H}/site*_block_ld_lmm_gen1.csv"))
    sites, cols_s, cols_r, nplots = [], {}, {}, {}
    base = None
    for f in files:
        sid = int(os.path.basename(f).split("_")[0][4:])
        d = pd.read_csv(f).set_index("gid")
        if base is None:
            base = d[["chrom", "unit_start", "unit_end", "panel_freq"]].copy()
        sites.append(sid); cols_s[sid] = d["s"]; cols_r[sid] = d["resid"]
        mj = json.load(open(f.replace(".csv", "_meta.json"))); nplots[sid] = mj["n_plot"]
    sites = sorted(sites)
    S = np.column_stack([cols_s[s].reindex(base.index).to_numpy() for s in sites])   # (B, nsite)
    R = np.column_stack([cols_r[s].reindex(base.index).to_numpy() for s in sites])
    w = np.array([nplots[s] for s in sites], float)
    print(f"{len(sites)} sites x {S.shape[0]:,} blocks; site weights (n_plot) {w.min():.0f}-{w.max():.0f}")

    # ---- climate (garden) ----
    cl = pd.read_csv(CLIM).set_index("site")
    bios = [c for c in cl.columns if c.startswith("bio")]
    cl = cl.loc[[s for s in sites if s in cl.index]]
    miss = [s for s in sites if s not in cl.index]
    if miss:
        print(f"  WARNING: {len(miss)} sites missing climate: {miss}")
    keep = [s for s in sites if s in cl.index]
    idx = [sites.index(s) for s in keep]
    S = S[:, idx]; R = R[:, idx]; w = w[idx]
    X = stats.zscore(cl[bios].to_numpy(), axis=0)
    pc1 = X @ np.linalg.svd(X, full_matrices=False)[2][0]          # climate PC1 (main gradient)
    bio1 = cl["bio1"].to_numpy()
    print(f"  {len(keep)} sites with climate; PC1 var-expl + bio1 (temp). "
          f"corr(PC1,bio1)={np.corrcoef(pc1,bio1)[0,1]:.2f}")

    En_S = wcorr_setup(S, w); En_R = wcorr_setup(R, w)
    out = base.copy()
    for cname, cvec in [("pc1", pc1), ("bio1", bio1)]:
        for elab, En in [("s", En_S), ("r", En_R)]:
            obs, pm, pf, q, thr = run_assoc(En, cvec, w, f"{elab} ~ {cname}")
            out[f"r_{elab}_{cname}"] = obs; out[f"pfwer_{elab}_{cname}"] = pf; out[f"q_{elab}_{cname}"] = q
    out.to_csv(f"{H}/crosssite_climate.csv", index=False)

    # ---- summary: ecotype-level (s) vs block-specific (r); Chr2:13.7Mb check ----
    print("\nSUMMARY  (climate-associated blocks, FWER p<0.05):")
    for cname in ["pc1", "bio1"]:
        ns = int((out[f"pfwer_s_{cname}"] < 0.05).sum()); nr = int((out[f"pfwer_r_{cname}"] < 0.05).sum())
        print(f"  {cname}:  s_raw (~ecotype) {ns}   |   r_resid (block-specific) {nr}")
    chr2 = out[(out.chrom == "Chr2") & (out.unit_start.between(13.6e6, 13.9e6))]
    if len(chr2):
        print("\n  Chr2:13.7Mb cluster — best climate assoc across its blocks:")
        for cn in ["pc1", "bio1"]:
            print(f"    {cn}: s |r|max={chr2[f'r_s_{cn}'].abs().max():.2f} (q_min={chr2[f'q_s_{cn}'].min():.3f}); "
                  f"r |r|max={chr2[f'r_r_{cn}'].abs().max():.2f} (q_min={chr2[f'q_r_{cn}'].min():.3f})")
    print(f"\n[done] {H}/crosssite_climate.csv")


if __name__ == "__main__":
    main()
