#!/usr/bin/env python
"""Latent-factor-adjusted PARTIAL RANK (Spearman) climate-GEA.

LFMM idea (K latent factors absorb population structure) applied to a monotone
rank statistic, fully vectorized across millions of variants.

Per variant i, the partial Spearman between the variant's AF and the climate
variable, adjusting out K structure axes U (LFMM's latent factors):

  * rank-transform env across pools -> re;  rank-transform each variant's AF -> Ra
  * design D = [1 | U]  (355 x (1+K)); projection  M_perp = I - D (D'D)^-1 D'
  * residualize:  re_resid = M_perp re ;  Ra_resid = Ra M_perp'
  * partial r_i = corr(Ra_resid[i], re_resid)
  * t_i = r_i sqrt(df/(1-r_i^2)),  df = n - K - 2 ; two-sided p from Student-t

K=0 recovers ordinary Spearman (unadjusted). Latent factors are the top-K
eigenvectors of the pool-pool gram matrix of the centered AF matrix, standardized
(identical to latent_factors() in run_binomial_latent.py).

Writes (--out):
  lf_rank_{class}_gen{g}_{climate}.csv   chrom,pos,ref_len,alt_len,MAF,block,r,pval,pval_unadj
  lf_rank_{class}_gen{g}_{climate}.npz   p (K-adj), p_unadj, r, gif curve, K, meta

Usage (kmate env):
  $PY run_lf_rank.py --class snp    --gen 9 --climate bio1 --k 16
  $PY run_lf_rank.py --class nonsnp --gen 9 --climate bio1 --k 16
"""
from __future__ import annotations
import argparse, os, sys, warnings
import numpy as np
import pandas as pd
from scipy.stats import t as tdist, chi2
from sklearn.preprocessing import StandardScaler
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

warnings.filterwarnings("ignore")


def rank_avg(A: np.ndarray) -> np.ndarray:
    """Average (tie-aware) ranks along axis=1, 1-based. A: [M, n] -> [M, n]."""
    A = np.atleast_2d(A)
    M, n = A.shape
    order = np.argsort(A, axis=1, kind="stable")
    sA = np.take_along_axis(A, order, axis=1)
    idx = np.arange(n)[None, :]
    # start index of each tie-group (in sorted order)
    is_new = np.ones((M, n), bool)
    is_new[:, 1:] = sA[:, 1:] != sA[:, :-1]
    start = np.where(is_new, idx, -1)
    np.maximum.accumulate(start, axis=1, out=start)
    # end index of each tie-group
    is_last = np.ones((M, n), bool)
    is_last[:, :-1] = sA[:, :-1] != sA[:, 1:]
    end = np.where(is_last, idx, n)
    end = np.minimum.accumulate(end[:, ::-1], axis=1)[:, ::-1]
    avg_sorted = (start + end) / 2.0 + 1.0                 # 1-based average rank
    ranks = np.empty((M, n), np.float64)
    np.put_along_axis(ranks, order, avg_sorted, axis=1)
    return ranks


def latent_factors(af: np.ndarray, k: int) -> np.ndarray:
    """Top-k per-pool structure loadings (LFMM's U) from [pools x rec] AF matrix."""
    Xc = af.astype(np.float64)
    Xc = Xc - Xc.mean(axis=0, keepdims=True)              # center each variant
    G = Xc @ Xc.T
    w, V = np.linalg.eigh(G)
    U = V[:, ::-1][:, :k]
    return StandardScaler().fit_transform(U)


def m_perp(U: np.ndarray | None, n: int) -> np.ndarray:
    """Residual projector I - D(D'D)^-1 D', D = [1 | U] (or just [1] if U is None)."""
    ones = np.ones((n, 1))
    D = ones if U is None else np.column_stack([ones, U])
    P = D @ np.linalg.pinv(D.T @ D) @ D.T
    return np.eye(n) - P


def partial_spearman(Ra: np.ndarray, Mp: np.ndarray, re_res: np.ndarray) -> np.ndarray:
    """Partial Spearman r for each row of ranks Ra vs residualized env re_res."""
    Ra_res = Ra @ Mp                                       # Mp symmetric
    num = Ra_res @ re_res
    den = np.linalg.norm(Ra_res, axis=1) * np.linalg.norm(re_res)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = num / den
    r[~np.isfinite(r)] = np.nan
    return r


def pval_from_r(r: np.ndarray, df: int) -> np.ndarray:
    rc = np.clip(r, -0.999999, 0.999999)
    tstat = rc * np.sqrt(df / (1.0 - rc * rc))
    p = 2.0 * tdist.sf(np.abs(tstat), df)
    p[~np.isfinite(r)] = np.nan
    return p


def gif(p: np.ndarray) -> float:
    p = p[np.isfinite(p)]
    if p.size == 0:
        return np.nan
    p = np.clip(p, 1e-300, 1.0)
    return float(np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", required=True,
                    choices=["snp", "sv", "smallindel", "nonsnp"])
    ap.add_argument("--gen", type=int, default=9)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--k", type=int, default=16)
    ap.add_argument("--cmdir", default=f"{lib.GEA}/phase1_replication/results/class_matrices")
    ap.add_argument("--out", default=f"{lib.GEA}/gea_newpanel/results/lf_rank")
    ap.add_argument("--chunk", type=int, default=100000)
    ap.add_argument("--kcurve", default="0,4,8,16,32")
    ap.add_argument("--curve-sub", type=int, default=200000,
                    help="subsample size for the GIF-vs-K curve")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    recs = pd.read_csv(f"{args.cmdir}/{args.cls}_gen{args.gen}.records.csv")
    pools = pd.read_csv(f"{args.cmdir}/gen{args.gen}.pools.csv")
    clim = pools[args.climate].to_numpy(float)
    n = len(clim)
    af = np.load(f"{args.cmdir}/{args.cls}_gen{args.gen}_af.npy")       # [pools x rec]
    assert af.shape[0] == n and af.shape[1] == len(recs)
    M = af.shape[1]
    K = args.k

    # latent factors (compute at max K needed; nested eigenvectors)
    kmax = max([K] + [int(x) for x in args.kcurve.split(",") if x])
    Umax = latent_factors(af, kmax) if kmax > 0 else np.zeros((n, 0))

    # env ranks, residualized under adjusted (K) and unadjusted (0) projectors
    re = rank_avg(clim)[0]
    Mp_adj = m_perp(Umax[:, :K] if K > 0 else None, n)
    Mp_un = m_perp(None, n)
    re_adj = Mp_adj @ re
    re_un = Mp_un @ re
    df_adj, df_un = n - K - 2, n - 2

    print(f"{args.cls} gen{args.gen}: {n} pools x {M:,} records vs {args.climate}; "
          f"K={K} (df_adj={df_adj}); ranking + partial Spearman", flush=True)

    r_adj = np.full(M, np.nan)
    p_adj = np.full(M, np.nan)
    p_un = np.full(M, np.nan)
    for a in range(0, M, args.chunk):
        b = min(a + args.chunk, M)
        Ra = rank_avg(af[:, a:b].T.astype(np.float64))                  # [chunk x n]
        r = partial_spearman(Ra, Mp_adj, re_adj)
        ru = partial_spearman(Ra, Mp_un, re_un)
        r_adj[a:b] = r
        p_adj[a:b] = pval_from_r(r, df_adj)
        p_un[a:b] = pval_from_r(ru, df_un)
        print(f"  {b:,}/{M:,}", end="\r", flush=True)
    print()

    # GIF-vs-K curve on a subsample (recompute ranks once for the subsample)
    ks = [int(x) for x in args.kcurve.split(",") if x]
    rng = np.random.default_rng(0)
    sub = rng.choice(M, size=min(args.curve_sub, M), replace=False)
    Ra_sub = rank_avg(af[:, sub].T.astype(np.float64))
    gif_curve = {}
    for k in ks:
        Mp = m_perp(Umax[:, :k] if k > 0 else None, n)
        rk = partial_spearman(Ra_sub, Mp, Mp @ re)
        gif_curve[k] = gif(pval_from_r(rk, n - k - 2))

    g_un, g_adj = gif(p_un), gif(p_adj)
    print(f"  GIF: unadj(K=0)={g_un:.3f}  adj(K={K})={g_adj:.3f}", flush=True)
    print(f"  GIF-vs-K: " + "  ".join(f"K{k}={gif_curve[k]:.2f}" for k in ks), flush=True)

    out_df = recs.rename(columns={"maf": "MAF"}).assign(r=r_adj, pval=p_adj, pval_unadj=p_un)
    cols = ["chrom", "pos", "ref_len", "alt_len", "MAF", "block", "r", "pval", "pval_unadj"]
    csv = f"{args.out}/lf_rank_{args.cls}_gen{args.gen}_{args.climate}.csv"
    out_df[cols].to_csv(csv, index=False)
    npz = f"{args.out}/lf_rank_{args.cls}_gen{args.gen}_{args.climate}.npz"
    np.savez_compressed(npz, p=p_adj, p_unadj=p_un, r=r_adj,
                        gif_adj=g_adj, gif_unadj=g_un, K=K,
                        gif_curve_k=np.array(ks),
                        gif_curve_v=np.array([gif_curve[k] for k in ks]),
                        chrom=recs["chrom"].to_numpy().astype("U8"),
                        pos=recs["pos"].to_numpy(), block=recs["block"].to_numpy().astype("U16"))

    fin = np.isfinite(p_adj)
    bonf = 0.05 / int(fin.sum())
    print(f"  tested {int(fin.sum()):,}/{M:,} | adj p<Bonf({bonf:.1e}): "
          f"{int((p_adj[fin] < bonf).sum()):,} | adj p<1e-5: {int((p_adj[fin] < 1e-5).sum()):,}"
          f" | min p={np.nanmin(p_adj):.2e}", flush=True)
    print(f"  -> {csv}\n  -> {npz}", flush=True)


if __name__ == "__main__":
    main()
