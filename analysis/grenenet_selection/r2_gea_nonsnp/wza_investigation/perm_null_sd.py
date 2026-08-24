#!/usr/bin/env python
"""Permutation-null test of the WZA null-SD-vs-SNP-count shape (isotonic justification).

Question: is the null SD of the block weighted-Z monotone-non-decreasing in block
SNP count (justifying an isotonic SD fit), or does it genuinely dip/plateau?

Null: SITE-permute the climate (bio1) across the 31 GrENE-Net sites (breaks the
genotype-climate association, preserves LD + pool/site structure), recompute a
per-SNP rank-correlation (Spearman — same rank family as the production Kendall
model; gen9 AF has 0 missing so this is exact & vectorized), then build the
PRODUCTION block statistic (genome-wide rank-transform of p -> z -> weighted-Z per
mcf90 block, cap-resampled) and its rolling SD-vs-n. Repeat NPERM times.

Also, under each null replicate, count how many blocks the deg-2 vs isotonic
correction sends to fabricated p==0 / NaN. deg-2 fabricating under a signal-free
null == proof the p==0 is a pure correction artifact, not real signal.

Saves perm_null_sd.npz for the notebook to plot. Run via sbatch (128G).
"""
import sys
import numpy as np, pandas as pd
from scipy.stats import norm, t as tdist
from sklearn.isotonic import IsotonicRegression

GEA = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection"
sys.path.insert(0, GEA)
import lib   # assign_clq_blocks (same mcf90 interval map reblock.py uses)
CM  = f"{GEA}/phase1_replication/results/class_matrices"
WIN = f"{GEA}/phase1_replication/results/clq90/wza_in"
CLASSES = ["snp", "smallindel", "sv"]
CAP = {"snp": 1000, "sv": 350, "smallindel": 350}
NPERM, MAF_MIN, ROLLER, MINE = 20, 0.05, 50, 40
NPOOL = 352

pools = pd.read_csv(f"{CM}/gen9.pools.csv")
site = pools["site"].to_numpy()
usites = np.unique(site); site_pos = {s: i for i, s in enumerate(usites)}
site_ix = np.array([site_pos[s] for s in site])
site_bio1 = pools.groupby("site")["bio1"].first().reindex(usites).to_numpy()  # per-site climate
print(f"{len(pools)} pools, {len(usites)} sites", flush=True)

def unitrank(v):  # rank, center, unit-norm a 1-D vector
    r = pd.Series(v).rank().to_numpy().astype(np.float64)   # copy (rank().to_numpy() is read-only)
    r = r - r.mean(); n = np.linalg.norm(r)
    return r / n if n else r

def block_stat(pval, pq, blocks, cap, seed):
    """Production block weighted-Z: genome-wide rank-transform p, z, cap-resample."""
    N = len(pval)
    ranks = pd.Series(pval).rank(method="first").to_numpy()
    pVal = np.clip(ranks / N, 1e-15, 1 - 1e-3)
    z = norm.ppf(1 - pVal)
    d = pd.DataFrame({"b": blocks, "num": pq * z, "den": pq ** 2})
    g = d.groupby("b"); Z = g["num"].sum() / np.sqrt(g["den"].sum()); cnt = g.size()
    res = pd.DataFrame({"SNPs_raw": cnt, "Z": Z})
    res = res[res.SNPs_raw >= 2].copy(); res["SNPs"] = res.SNPs_raw.clip(upper=cap)
    big = res.index[res.SNPs_raw > cap]
    if len(big):
        rng = np.random.default_rng(seed)
        bser = pd.Series(np.arange(len(blocks))).groupby(blocks).apply(lambda s: s.to_numpy())
        for b in big:
            ix = bser[b]; pqb = pq[ix]; zb = z[ix]; n = len(ix)
            res.loc[b, "Z"] = np.mean([(pqb[c] * zb[c]).sum() / np.sqrt((pqb[c] ** 2).sum())
                                       for c in (rng.choice(n, cap, replace=False) for _ in range(100))])
    return res

def roll_sd(res):
    s = res[res.Z.notna()].sort_values("SNPs")
    var = s.Z.rolling(ROLLER, min_periods=MINE).var()
    x = s.SNPs.rolling(ROLLER, min_periods=MINE).mean()
    m = var.notna()
    return x[m].to_numpy(), np.sqrt(var[m].to_numpy())

def fab_counts(res, x, sd):
    """#p==0 and #NaN under deg2 vs isotonic, for a given block set + rolling SD."""
    n = res.SNPs.to_numpy(); out = {}
    for how in ("deg2", "isotonic"):
        if how == "deg2":
            sdp = np.poly1d(np.polyfit(x, sd, 2))(n)
        else:
            ir = IsotonicRegression(increasing=True, out_of_bounds="clip"); ir.fit(x, sd); sdp = ir.predict(n)
        with np.errstate(invalid="ignore"):
            p = 1 - norm.cdf(res.Z.to_numpy(), loc=0.0, scale=sdp)  # null mean ~0
        out[how] = (int((p == 0).sum()), int(np.isnan(p).sum()))
    return out

GRID = {c: np.linspace(2, CAP[c], 120) for c in CLASSES}
save = {}
for cls in CLASSES:
    print(f"\n=== {cls} ===", flush=True)
    af = np.load(f"{CM}/{cls}_gen9_af.npy")                      # [pools x recs], row-aligned to records.csv (no NaN)
    rec = pd.read_csv(f"{CM}/{cls}_gen9.records.csv")            # chrom,pos,maf ; af column = record ROW index
    # mcf90 block straight from position (same interval map + r2 reblock.py uses); drop gaps + MAF filter
    rec["block"] = lib.assign_clq_blocks(rec["chrom"].to_numpy(str), rec["pos"].to_numpy(np.int64), r2=0.9)
    rec["afrow"] = np.arange(len(rec))
    keep = (rec["block"].to_numpy() != "") & (rec["maf"].to_numpy() > MAF_MIN)
    sub = rec[keep]
    cols = sub["afrow"].to_numpy(); blocks = sub["block"].to_numpy().astype(str)
    pq = (sub["maf"] * (1 - sub["maf"])).to_numpy()
    afx = af[:, cols].astype(np.float64)                        # [pools x inblock]
    print(f"  in-block MAF>0.05 records: {len(sub):,}", flush=True)
    # rank AF across pools once (Spearman ranks), center + unit-norm per column
    rk = afx.argsort(0).argsort(0).astype(np.float64) + 1.0
    rk -= rk.mean(0); nrm = np.sqrt((rk ** 2).sum(0)); nrm[nrm == 0] = 1.0; rk /= nrm
    del afx
    n = NPOOL

    def spearman_p(clim_pool):
        xr = unitrank(clim_pool)                                # [pools]
        r = xr @ rk                                             # [recs] spearman
        r = np.clip(r, -0.999999, 0.999999)
        tt = r * np.sqrt((n - 2) / (1 - r ** 2))
        return 2 * tdist.sf(np.abs(tt), n - 2)

    # REAL (observed climate) curve for overlay
    p_real = spearman_p(pools["bio1"].to_numpy())
    res_r = block_stat(p_real, pq, blocks, CAP[cls], seed=0)
    xr_, sdr = roll_sd(res_r); real_curve = np.interp(GRID[cls], xr_, sdr, left=np.nan, right=sdr[-1])

    null_curves = np.full((NPERM, len(GRID[cls])), np.nan)
    fab = []
    for r_ in range(NPERM):
        rng = np.random.default_rng(1000 + r_)
        perm_site_bio1 = site_bio1[rng.permutation(len(usites))]  # shuffle climate across sites
        clim = perm_site_bio1[site_ix]                            # broadcast to pools
        pv = spearman_p(clim)
        res = block_stat(pv, pq, blocks, CAP[cls], seed=r_)
        xx, ss = roll_sd(res)
        null_curves[r_] = np.interp(GRID[cls], xx, ss, left=np.nan, right=ss[-1])
        fc = fab_counts(res, xx, ss)
        fab.append((fc["deg2"][0], fc["deg2"][1], fc["isotonic"][0], fc["isotonic"][1]))
    fab = np.array(fab)
    # monotonicity of the mean null curve up to the support edge
    mean_null = np.nanmean(null_curves, 0)
    ok = np.isfinite(mean_null)
    diffs = np.diff(mean_null[ok]); frac_up = float((diffs >= -1e-9).mean())
    print(f"  NULL deg2 fabricated p==0/rep: mean={fab[:,0].mean():.1f} max={fab[:,0].max()}  "
          f"NaN/rep mean={fab[:,1].mean():.1f}", flush=True)
    print(f"  NULL isotonic p==0/rep: mean={fab[:,2].mean():.1f}  NaN/rep mean={fab[:,3].mean():.1f}", flush=True)
    print(f"  mean null SD monotone-nondecreasing fraction (up support): {frac_up:.3f}", flush=True)
    save[f"{cls}_grid"] = GRID[cls]; save[f"{cls}_real"] = real_curve
    save[f"{cls}_null"] = null_curves; save[f"{cls}_fab"] = fab

np.savez(f"{GEA}/wza_investigation/perm_null_sd.npz", **save)
print("\nwrote perm_null_sd.npz", flush=True)
