#!/usr/bin/env python
"""K sweep for the latent-factor-adjusted partial-rank (Spearman) climate-GEA:
find the K that de-inflates (GIF->1) WITHOUT killing signal. For each K and class
(snp,nonsnp) and a few axes (bio1, bio18, pc1), report GIF and hit counts under the
K-adjusted p (raw) and after GIF-calibration. Writes a json + a GIF-vs-K / hits-vs-K
figure so we can pick the knee instead of defaulting to K=16.

Efficiency: rank each variant's AF ONCE; latent factors nested to Kmax; per K the
residual projector Mp changes, so Ra_res = Ra@Mp is computed once per K and reused
across axes. SNP-derived factors (common structure reference).
Run in kmate env; plot separately in basic.
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
from scipy.stats import chi2
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.dirname(HERE))
import lib
from run_lf_rank import rank_avg, latent_factors, m_perp, pval_from_r, gif

CM = f"{lib.GEA}/phase1_replication/results/class_matrices"
OUT = f"{lib.GEA}/gea_newpanel/results/k_sweep"
KGRID = [0, 2, 4, 6, 8, 10, 12, 16, 20, 24, 32]
AXES = ["bio1", "bio18", "pc1"]
GEN = 9


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = (p[o] * n) / (np.arange(n) + 1); q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.clip(q, 0, 1)


def calibrate(p, g):
    if not np.isfinite(g) or g < 1.0:
        return p
    return chi2.sf(chi2.isf(np.clip(p, 1e-300, 1.0), 1) / g, 1)


def main():
    os.makedirs(OUT, exist_ok=True)
    pools = pd.read_csv(f"{CM}/gen{GEN}.pools.csv")
    n = len(pools)
    bios = [f"bio{i}" for i in range(1, 20)]
    Bz = pools[bios].to_numpy(float)
    Bz = (Bz - Bz.mean(0)) / Bz.std(0, ddof=0)
    _, _, Vt = np.linalg.svd(Bz - Bz.mean(0), full_matrices=False)
    env = {ax: pools[ax].to_numpy(float) for ax in bios}
    env["pc1"] = Bz @ Vt[0]
    re = {ax: rank_avg(env[ax])[0] for ax in AXES}

    kmax = max(KGRID)
    snp_af = np.load(f"{CM}/snp_gen{GEN}_af.npy")
    U = latent_factors(snp_af, kmax)
    del snp_af

    rows = []
    for cls in ["snp", "nonsnp"]:
        recs = pd.read_csv(f"{CM}/{cls}_gen{GEN}.records.csv")
        af = np.load(f"{CM}/{cls}_gen{GEN}_af.npy")
        M = af.shape[1]
        Ra = np.empty((M, n), np.float32)
        for a in range(0, M, 100000):
            b = min(a + 100000, M)
            Ra[a:b] = rank_avg(af[:, a:b].T.astype(np.float64)).astype(np.float32)
        del af
        for K in KGRID:
            Mp = m_perp(U[:, :K] if K > 0 else None, n).astype(np.float32)
            Ra_res = Ra @ Mp                                   # [M x n]
            ra_norm = np.linalg.norm(Ra_res, axis=1); ra_norm[ra_norm == 0] = np.inf
            df = n - K - 2
            for ax in AXES:
                rr = (Mp @ re[ax].astype(np.float32))
                r = (Ra_res @ rr) / (ra_norm * np.linalg.norm(rr))
                p = pval_from_r(r.astype(np.float64), df)
                g = gif(p)
                bonf = 0.05 / M
                pc = calibrate(p, g)
                rows.append(dict(cls=cls, K=K, axis=ax, GIF=round(float(g), 3),
                    min_p=float(np.nanmin(p)),
                    n_bonf_raw=int((p < bonf).sum()),
                    n_p1e5_raw=int((p < 1e-5).sum()),
                    n_bonf_gifcal=int((pc < bonf).sum()),
                    n_fdr_gifcal=int((bh(pc) < 0.05).sum())))
            del Ra_res
            print(f"[{cls} K={K:>2}] " + " | ".join(
                f"{ax} GIF={[x for x in rows if x['cls']==cls and x['K']==K and x['axis']==ax][0]['GIF']:.2f}"
                f" bonf_raw={[x for x in rows if x['cls']==cls and x['K']==K and x['axis']==ax][0]['n_bonf_raw']}"
                for ax in AXES), flush=True)
        del Ra
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/k_sweep.csv", index=False)
    json.dump(rows, open(f"{OUT}/k_sweep.json", "w"), indent=2)
    print(f"\nwrote {OUT}/k_sweep.csv")


if __name__ == "__main__":
    main()
