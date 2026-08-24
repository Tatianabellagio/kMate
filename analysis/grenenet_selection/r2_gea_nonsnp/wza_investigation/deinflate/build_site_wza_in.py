#!/usr/bin/env python
"""Build SITE-level LFMM per-record p as a drop-in replacement for the pool-level
wza_in_clq09_tile files, so the raw-manhattan notebook can plot the site-collapsed scan.

Site collapse = mean Δp across the ~11 pools within each of the 31 sites (removes the
pseudoreplication that inflates the pool-level scan), then lfmm_ridge(K=1) + lfmm_test.
K=1 is the site-level operating point: at n=31 sites K=1 minimized the GIF (higher K
re-inflates by overfitting the 31 samples).

Writes results/multiaxis/wza_in_clq09_tile_site/lfmm_{cls}_gen9_{axis}.csv with the SAME
schema/rows as the pool file (chrom,pos,ref_len,alt_len,MAF,block) but:
    pval = site-level RAW lfmm p   (lfmm_test calibrate=NULL -- NO genomic-inflation correction)
Usage: build_site_wza_in.py [--classes snp,sv] [--K 1]
"""
from __future__ import annotations
import argparse, os, subprocess, sys, tempfile
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
import lib

HERE = os.path.dirname(os.path.abspath(__file__))
CM = f"{lib.GEA}/phase1_replication/results/class_matrices"
POOL_WZAIN = f"{lib.GEA}/phase1_replication/results/multiaxis/wza_in_clq09_tile"
OUT = f"{lib.GEA}/phase1_replication/results/multiaxis/wza_in_clq09_tile_site"
RSCRIPT = "/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript"
RUNNER = f"{HERE}/run_lfmm_scores_multiaxis.R"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
_P0 = {"snp": "p0_snp.npy", "sv": "p0_nonsnp.npy", "smallindel": "p0_nonsnp.npy", "nonsnp": "p0_nonsnp.npy"}


def build_site_Y(cls, stem):
    """Site-mean Δp matrix [n_site x n_rec] + dims; returns (recs, usites)."""
    recs = pd.read_csv(f"{CM}/{cls}_gen9.records.csv")
    pools = pd.read_csv(f"{CM}/gen9.pools.csv")
    af = np.load(f"{CM}/{cls}_gen9_af.npy")
    p0 = np.load(f"{lib.AF_STORE}/{_P0[cls]}")[recs["col"].to_numpy()]
    dp = af - p0[None, :]
    with np.errstate(invalid="ignore"):
        colmean = np.nanmean(dp, axis=0)
    colmean = np.where(np.isfinite(colmean), colmean, 0.0)
    miss = ~np.isfinite(dp); dp[miss] = np.take(colmean, np.where(miss)[1])
    site = pools["site"].to_numpy(); usites = pd.unique(site)
    S = np.vstack([dp[site == s].mean(0) for s in usites])       # [n_site x n_rec]
    S.astype(np.float64).tofile(f"{stem}_Y.f64")
    open(f"{stem}_dims.txt", "w").write(f"{len(usites)} {S.shape[1]}\n")
    for ax in AXES:
        xc = pools.groupby("site")[ax].first().reindex(usites).to_numpy(float)
        z = (xc - xc.mean()) / xc.std(ddof=0)
        pd.DataFrame({ax: z}).to_csv(f"{stem}_env_{ax}.csv", index=False)
    print(f"[{cls}] site Y {S.shape} ({len(pools)} pools -> {len(usites)} sites), env written for {len(AXES)} axes",
          flush=True)
    return recs, usites


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", default="snp,sv")
    ap.add_argument("--K", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    indir = f"{HERE}/inputs"; os.makedirs(indir, exist_ok=True)

    for cls in args.classes.split(","):
        stem = f"{indir}/lfmm_{cls}_site_gen9"
        recs, _ = build_site_Y(cls, stem)
        keys = ["chrom", "pos", "ref_len", "alt_len"]
        rk = recs[keys].reset_index(drop=True)

        with tempfile.TemporaryDirectory() as td:
            subprocess.run([RSCRIPT, RUNNER, stem, str(args.K), td, ",".join(AXES)], check=True)
            for ax in AXES:
                sc = pd.read_csv(f"{td}/scores_{ax}.csv")
                assert len(sc) == len(rk), (cls, ax, len(sc), len(rk))
                # pool wza_in is in class-matrix record order (== recs order == scores order);
                # variant keys are NOT unique (multiallelic), so assign POSITIONALLY after
                # asserting exact row alignment of the keys.
                pool = pd.read_csv(f"{POOL_WZAIN}/lfmm_{cls}_gen9_{ax}.csv")
                assert len(pool) == len(rk), (cls, ax, len(pool), len(rk))
                assert (pool[keys].to_numpy() == rk[keys].to_numpy()).all(), (cls, ax, "row misalignment")
                out = pool.copy()
                out["pval"] = sc["praw"].to_numpy()          # site-level RAW p (lfmm calibrate=NULL, no gif)
                out.to_csv(f"{OUT}/lfmm_{cls}_gen9_{ax}.csv", index=False)
        print(f"[{cls}] wrote {len(AXES)} site wza_in files -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
