#!/usr/bin/env python
"""Per-block fine-map inputs for the driver-vs-passenger test (SuSiE-RSS).

For each target clq0.9 block (default: the SV-containing BH-sig blocks from the
all-class WZA screen), build the two SuSiE-RSS inputs on the SAME allele coding:

  * z  : SITE-LEVEL marginal association z per variant (all classes). Plots are
         flower-weighted-averaged to their site (N=31 sites = the effective unit;
         this is where the site-level effective-N enters), then z = t-stat of the
         Pearson correlation between site AF and site bio1. Using the site as the
         unit deflates the 355-pool pseudo-replication (GIF ~2.5) at the source.
  * R  : founder-panel signed LD (correlation of 231-founder alt dosage) among the
         block's variants. Full-rank reference LD (a 31-site in-sample LD would be
         rank<=31 -> degenerate for big blocks). Same alt allele as z -> signs match.

Driver readout downstream: is any SV in a 95% credible set / the top-PIP variant,
or is the SV PIP always <= the best SNP? That's the LD-confound-beating question.

Outputs (--out, default analysis/grenenet_selection/extras/driver_passenger/results/finemap_inputs):
  <block>.npz  keys: block, chrom, pos, ref_len, alt_len, cls, is_sv, beta, se, z,
                     R (m x m float32), sites (31,), matched (bool per variant)

Usage:
  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY build_finemap_inputs.py            # 22 SV-containing BH-sig blocks
  $PY build_finemap_inputs.py --blocks all_bh   # all 71 BH-sig blocks
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
from scipy import sparse

HERE = os.path.dirname(os.path.abspath(__file__))
GEA_DIR = os.path.dirname(HERE)
sys.path.insert(0, GEA_DIR)
import lib

CLASSES = ("snp", "smallindel", "sv")
CM = f"{lib.GEA}/phase1_replication/results/class_matrices"
PANEL = "/global/scratch/users/tbellg/kmate/panel/arch3"
DP = f"{lib.GEA}/driver_passenger/results"


def site_weights(pools: pd.DataFrame):
    """31 x 355 flower-weighted plot->site averaging matrix + per-site bio1."""
    sites = np.sort(pools["site"].unique())
    W = np.zeros((len(sites), len(pools)), dtype=np.float64)
    for i, s in enumerate(sites):
        m = (pools["site"].to_numpy() == s)
        w = pools["total_flowers"].to_numpy()[m].astype(float)
        w = w / w.sum() if w.sum() > 0 else np.full(m.sum(), 1.0 / m.sum())
        W[i, np.where(m)[0]] = w
    bio1 = np.array([pools.loc[pools["site"] == s, "bio1"].iloc[0] for s in sites])
    return sites, W, bio1


def site_level_z(site_af: np.ndarray, bio1: np.ndarray):
    """Per-variant t-stat (=beta/se) of Pearson(site_af, bio1) across sites."""
    n = len(bio1)
    x = bio1 - bio1.mean()
    sxx = (x * x).sum()
    Y = site_af - site_af.mean(axis=0, keepdims=True)          # sites x m
    sxy = x @ Y                                                # m
    beta = sxy / sxx
    resid = Y - np.outer(x, beta)
    sse = (resid ** 2).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        se = np.sqrt(sse / (n - 2) / sxx)
        z = beta / se
    z[~np.isfinite(z)] = 0.0                                    # invariant/degenerate -> null
    return beta, se, z


def load_class_records(target_blocks: set):
    """All variants (any class) in the target blocks, with af-column index."""
    rows = []
    for cls in CLASSES:
        r = pd.read_csv(f"{CM}/{cls}_gen9.records.csv")
        r["afcol"] = np.arange(len(r))                          # af.npy col == records row
        blk = lib.assign_clq_blocks(np.asarray(r["chrom"], dtype=str),
                                    np.asarray(r["pos"]), r2=0.9)
        r["block"] = blk
        r = r[r["block"].isin(target_blocks)].copy()
        r["cls"] = cls
        rows.append(r[["chrom", "pos", "ref_len", "alt_len", "cls", "afcol", "block"]])
        print(f"  {cls:11s}: {len(r):,} records in target blocks", flush=True)
    return pd.concat(rows, ignore_index=True)


def founder_matrix(chrom: str):
    """231 x R alt dosage (uncalled -> col mean) + (pos,ref_len,alt_len)->col lookup."""
    ci = chrom.replace("Chr", "")
    d = f"{PANEL}/chr{ci}"
    meta = dict(np.load(f"{d}/var_pa_231_arch3_chr{ci}.meta.npz", allow_pickle=True))
    vp = np.load(f"{d}/var_pa_231_arch3_chr{ci}.var_pa.npz", allow_pickle=True)
    vc = np.load(f"{d}/var_pa_231_arch3_chr{ci}.var_called.npz", allow_pickle=True)
    A = sparse.csr_matrix((vp["data"], vp["indices"], vp["indptr"]), shape=tuple(vp["shape"]))
    C = sparse.csr_matrix((vc["data"], vc["indices"], vc["indptr"]), shape=tuple(vc["shape"]))
    # orient founders x records, keep sparse (slice columns per block, densify small)
    if A.shape[0] != 231:
        A, C = A.T, C.T
    A, C = A.tocsc(), C.tocsc()
    key = {}
    for j, (p, rl, al) in enumerate(zip(meta["pos"], meta["ref_len"], meta["alt_len"])):
        key.setdefault((int(p), int(rl), int(al)), j)           # first wins on collision
    return A, C, key


def block_R(A, C, key, recs: pd.DataFrame):
    """Signed founder LD (m x m) for a block's variants; unmatched -> R row/col 0, matched flag."""
    m = len(recs)
    cols = np.full(m, -1, dtype=np.int64)
    for i, (_, r) in enumerate(recs.iterrows()):
        cols[i] = key.get((int(r["pos"]), int(r["ref_len"]), int(r["alt_len"])), -1)
    matched = cols >= 0
    R = np.eye(m, dtype=np.float32)
    if matched.sum() >= 2:
        idx = cols[matched]
        G = np.asarray(A[:, idx].todense(), dtype=np.float32)   # 231 x k
        M = np.asarray(C[:, idx].todense()).astype(bool)
        # impute uncalled to per-variant called mean
        for j in range(G.shape[1]):
            if (~M[:, j]).any():
                cm = G[M[:, j], j].mean() if M[:, j].any() else 0.0
                G[~M[:, j], j] = cm
        Gc = G - G.mean(axis=0, keepdims=True)
        sd = Gc.std(axis=0)
        good = sd > 0
        Rk = np.eye(matched.sum(), dtype=np.float32)
        if good.sum() >= 2:
            Cc = np.corrcoef(Gc[:, good].T).astype(np.float32)
            gi = np.where(good)[0]
            Rk[np.ix_(gi, gi)] = np.nan_to_num(Cc, nan=0.0)
        mi = np.where(matched)[0]
        R[np.ix_(mi, mi)] = Rk
    return R, matched


def af_corr(af_cols: np.ndarray):
    """Signed correlation (m x m) among variants' pool-AF vectors; invariant -> diag."""
    m = af_cols.shape[1]
    R = np.eye(m, dtype=np.float32)
    sd = af_cols.std(axis=0)
    good = sd > 0
    if good.sum() >= 2:
        C = np.corrcoef(af_cols[:, good].T).astype(np.float32)
        gi = np.where(good)[0]
        R[np.ix_(gi, gi)] = np.nan_to_num(C, nan=0.0)
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", default="sv_bh",
                    help="'sv_bh' (22 SV-containing BH-sig, default) | 'all_bh' (71) | path to csv w/ block col")
    ap.add_argument("--level", default="site", choices=["site", "pool"],
                    help="site: 31-site z + founder LD (honest primary); "
                         "pool: 355-pool z + pool-AF LD (consistent full-power sensitivity)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if args.out is None:
        args.out = f"{DP}/finemap_inputs" + ("_pool" if args.level == "pool" else "")
    os.makedirs(args.out, exist_ok=True)

    comp = pd.read_csv(f"{DP}/block_composition_kendall_gen9_bio1.csv")
    comp = comp[comp["Z_pVal"].notna()].copy()
    comp["q"] = lib.bh(comp["Z_pVal"].to_numpy())
    if args.blocks == "sv_bh":
        tgt = comp[(comp["q"] < 0.05) & (comp["n_sv"] > 0)]["block"].tolist()
    elif args.blocks == "all_bh":
        tgt = comp[comp["q"] < 0.05]["block"].tolist()
    else:
        tgt = pd.read_csv(args.blocks)["block"].tolist()
    target_blocks = set(tgt)
    print(f"[fine-map inputs] {len(target_blocks)} target blocks", flush=True)

    pools = pd.read_csv(f"{CM}/gen9.pools.csv")
    sites, W, bio1 = site_weights(pools)
    pool_bio1 = pools["bio1"].to_numpy()
    print(f"  level={args.level}: "
          + ("31-site z (flower-wgt) + founder LD" if args.level == "site"
             else f"{len(pools)}-pool z + pool-AF LD (consistent, n={len(pools)})"), flush=True)

    recs = load_class_records(target_blocks)
    afmm = {c: np.load(f"{CM}/{c}_gen9_af.npy", mmap_mode="r") for c in CLASSES}

    done = 0
    for chrom, cg in recs.groupby("chrom"):
        A = C = key = None
        if args.level == "site":
            A, C, key = founder_matrix(chrom)
        for block, bg in cg.groupby("block"):
            bg = bg.sort_values("pos").reset_index(drop=True)
            # site-AF matrix (31 x m) via flower-weighted aggregation
            af_cols = np.empty((len(pools), len(bg)), dtype=np.float32)
            for c in CLASSES:
                sel = bg["cls"] == c
                if sel.any():
                    af_cols[:, np.where(sel.to_numpy())[0]] = afmm[c][:, bg.loc[sel, "afcol"].to_numpy()]
            if args.level == "site":
                obs = W @ af_cols                                # 31 x m (flower-wgt agg)
                yy, nobs = bio1, len(sites)
                beta, se, z = site_level_z(obs, yy)
                R, matched = block_R(A, C, key, bg)              # founder-panel LD
            else:                                                # pool: 355 obs, consistent LD
                obs = af_cols                                    # 355 x m
                yy, nobs = pool_bio1, len(pools)
                beta, se, z = site_level_z(obs, yy)
                R = af_corr(af_cols)                             # pool-AF LD (same data as z)
                matched = np.ones(len(bg), dtype=bool)
            is_sv = (bg["cls"] == "sv").to_numpy()
            np.savez_compressed(
                f"{args.out}/{block}.npz",
                block=block, chrom=chrom, pos=bg["pos"].to_numpy(),
                ref_len=bg["ref_len"].to_numpy(), alt_len=bg["alt_len"].to_numpy(),
                cls=bg["cls"].to_numpy(), is_sv=is_sv,
                beta=beta, se=se, z=z, R=R, matched=matched, sites=sites, n_obs=nobs)
            done += 1
            print(f"  {block:12s} m={len(bg):>4} (sv={int(is_sv.sum())}) "
                  f"matched={int(matched.sum())}/{len(bg)} max|z|={np.abs(z).max():.1f} "
                  f"sv_max|z|={np.abs(z[is_sv]).max() if is_sv.any() else 0:.1f}", flush=True)
    print(f"\n[done] {done} blocks -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
