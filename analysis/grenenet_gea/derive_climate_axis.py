#!/usr/bin/env python
"""Derive a climate-axis variant of a finished multisite founder-GWAS run WITHOUT re-running EMMAX.

The Bolormaa CLIMATE contrast is a pure function of the saved (Z, C, climate-vector): JOINT and GLOBAL
are climate-independent, so only the CLIMATE column, its permutation null, and the clade-level
correlation change when we swap the climate axis (e.g. bio1 -> bioclim PC1). This reads the finished
`--in-suffix` artifacts and writes parallel `--out-suffix` artifacts with the CLIMATE part recomputed on
the chosen axis, so the standard notebook builder renders an identical notebook keyed to that axis.

Writes (for OUT): multisite_founder_gwas{OUT}.{csv,npz,_meta.json}, multisite_climate_perm{OUT}.json,
cross_site_winners_30{OUT}.{csv,json}. Env: kmate.  Usage:
  python derive_climate_axis.py --axis bioPC1 --in-suffix _clq90 --out-suffix _clq90_pc1
"""
import os, sys, json, argparse
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

H = "results/grenenet_gea/hapfreq"
N_PERM, SEED, THR = 10000, 0, 2.5


def bh(pv):
    m = len(pv); o = pv.argsort(); q = np.empty(m)
    q[o] = np.minimum.accumulate((pv[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def climate_vector(axis, sites, bio1):
    """Return (per-site climate scores, human label) for the requested axis."""
    if axis == "bio1":
        return bio1.astype(float), "bio1 (mean annual temperature)"
    clim = lib.load_climate().reindex(sites)
    if axis == "bio12":
        return clim["bio12"].to_numpy(float), "bio12 (annual precipitation)"
    if axis == "bioPC1":
        B = clim[[f"bio{i}" for i in range(1, 20)]].to_numpy(float)
        Bz = (B - B.mean(0)) / B.std(0)
        U, sv, _ = np.linalg.svd(Bz - Bz.mean(0), full_matrices=False)
        pc1 = U[:, 0] * sv[0]
        if np.corrcoef(pc1, bio1)[0, 1] < 0:
            pc1 = -pc1
        ve = sv[0] ** 2 / (sv ** 2).sum()
        return pc1, f"bioclim PC1 ({ve:.0%} var, temp-precip composite)"
    raise SystemExit(f"unknown axis {axis}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--axis", default="bioPC1")
    ap.add_argument("--in-suffix", default="_clq90")
    ap.add_argument("--out-suffix", default="_clq90_pc1")
    a = ap.parse_args()
    I, O = a.in_suffix, a.out_suffix

    d = np.load(f"{H}/multisite_founder_gwas{I}.npz", allow_pickle=True)
    Z, C, sites, bio1 = d["Z"], d["C"], d["sites"], d["bio1"]
    M, S = Z.shape
    meta = json.load(open(f"{H}/multisite_founder_gwas{I}_meta.json"))
    cvar, label = climate_vector(a.axis, sites, bio1)
    print(f"[derive] axis={a.axis} ({label}) | {M:,} blocks x {S} sites | in={I} out={O}")

    # ---- CLIMATE contrast on the new axis (C-orthogonalized against GLOBAL 1-vector) ----
    Cinv = np.linalg.pinv(C); one = np.ones(S); dg = float(one @ Cinv @ one); ZC = Z @ Cinv

    def clim_z(b):
        c0 = (b - b.mean()) / b.std(); c = c0 - (float(one @ Cinv @ c0) / dg) * one
        return (ZC @ c) / np.sqrt(float(c @ Cinv @ c))

    z = clim_z(cvar); p = 2 * stats.norm.sf(np.abs(z)); q = bh(p)
    lam = float(np.median(z ** 2) / stats.chi2.ppf(0.5, 1)); bonf = 0.05 / M
    # Everything is rebuilt in npz-row (coordinate) order from Z/C + the coordinate arrays, so no
    # dependence on the in-csv row order. JOINT/GLOBAL are climate-independent -> recomputed identically.
    chrom = d["chrom"].astype(str); start = d["start"]; end = d["end"]; mac = d["mac"]
    Q = np.einsum("mi,ij,mj->m", Z, Cinv, Z); p_joint = stats.chi2.sf(Q, S); q_joint = bh(p_joint)
    ggl = (Z @ Cinv @ one) / np.sqrt(dg); p_global = 2 * stats.norm.sf(np.abs(ggl)); q_global = bh(p_global)

    out = pd.DataFrame({"chrom": chrom, "start": start, "end": end,
                        "unit": [f"{c}:{s}-{e}" for c, s, e in zip(chrom, start, end)], "mac": mac,
                        "chi2_joint": Q, "p_joint": p_joint, "q_joint": q_joint,
                        "z_global": ggl, "p_global": p_global, "q_global": q_global,
                        "z_clim": z, "p_clim": p, "q_clim": q})
    out.sort_values("p_joint").to_csv(f"{H}/multisite_founder_gwas{O}.csv", index=False)
    np.savez(f"{H}/multisite_founder_gwas{O}.npz", Z=Z, sites=sites, bio1=cvar, C=C,
             chrom=chrom, start=start, end=end, mac=mac)

    com = mac >= meta["mac_grm"]
    meta_o = dict(meta)
    meta_o["climate_axis"] = a.axis; meta_o["climate_label"] = label
    meta_o["tests"] = dict(meta["tests"])
    meta_o["tests"]["CLIMATE"] = dict(lam=lam, n_bonf=int((p < bonf).sum()), n_fdr=int((q < 0.05).sum()),
                                      n_fdr_common=int((q < 0.05)[com].sum()), best_q=float(q.min()))
    json.dump(meta_o, open(f"{H}/multisite_founder_gwas{O}_meta.json", "w"), indent=2)

    # ---- permutation null on the new axis (FWER max-stat + diffuse excess) ----
    obs = z; obs_max = float(np.abs(obs).max()); obs_nexc = int((np.abs(obs) > THR).sum()); obs_var = float(obs.var())
    rng = np.random.default_rng(SEED)
    nmax = np.empty(N_PERM); nexc = np.empty(N_PERM, int); nvar = np.empty(N_PERM)
    for k in range(N_PERM):
        zk = clim_z(rng.permutation(cvar)); ak = np.abs(zk)
        nmax[k] = ak.max(); nexc[k] = int((ak > THR).sum()); nvar[k] = zk.var()
    # clade winning-clade PC1 vs the new axis (independent 30-site file; win_pc1 is axis-independent)
    csw = pd.read_csv(f"{H}/cross_site_winners_30{I}.csv")
    cj_in = json.load(open(f"{H}/cross_site_winners_30{I}.json"))
    win = csw.set_index("site").reindex(sites)["win_pc1"].to_numpy(float)
    fin = np.isfinite(cvar) & np.isfinite(win)
    r = float(np.corrcoef(cvar[fin], win[fin])[0, 1])
    rng2 = np.random.default_rng(SEED)
    nd = np.array([abs(np.corrcoef(rng2.permutation(cvar[fin]), win[fin])[0, 1]) for _ in range(N_PERM)])
    clade_p = float((nd >= abs(r)).mean())

    perm = dict(n_perm=N_PERM, seed=SEED, n_blocks=int(M), n_sites=int(S), thr=THR,
                obs_max_z=obs_max, fwer_thr_95=float(np.quantile(nmax, 0.95)),
                n_fwer_sig=int((np.abs(obs) > np.quantile(nmax, 0.95)).sum()),
                fwer_p_top=float((nmax >= obs_max).mean()),
                obs_n_excess=obs_nexc, null_n_excess_med=float(np.median(nexc)),
                p_excess=float((nexc >= obs_nexc).mean()),
                obs_var_z=obs_var, null_var_z_med=float(np.median(nvar)), p_var=float((nvar >= obs_var).mean()),
                crosssite_pc1_bio1_r=r, crosssite_pc1_bio1_perm_p=clade_p,
                climate_axis=a.axis, note=f"CLIMATE axis = {label}; derived from {I} without re-running EMMAX")
    json.dump(perm, open(f"{H}/multisite_climate_perm{O}.json", "w"), indent=2)

    # clade artifacts (bio1 column carries the chosen axis so the notebook scatter x = this axis)
    csw_o = csw.copy(); csw_o["bio1"] = csw.site.map(dict(zip(sites.tolist(), cvar)))
    csw_o.to_csv(f"{H}/cross_site_winners_30{O}.csv", index=False)
    json.dump(dict(n_sites=int(fin.sum()), n_perm=N_PERM, mean_crosssite_r=cj_in.get("mean_crosssite_r"),
                   trait_gens=cj_in.get("trait_gens"), excluded=cj_in.get("excluded"),
                   pc1=dict(r=r, perm_p=clade_p), climate_axis=a.axis),
              open(f"{H}/cross_site_winners_30{O}.json", "w"), indent=2)

    print(f"[derive] CLIMATE: lam={lam:.2f} best_q={q.min():.3g} nFDR={int((q<0.05).sum())} | "
          f"FWER p(top)={perm['fwer_p_top']:.3f} excess p={perm['p_excess']:.3f} | "
          f"clade r={r:.2f} perm p={clade_p:.3f}")
    print(f"[done] wrote *{O}.* artifacts")


if __name__ == "__main__":
    main()
