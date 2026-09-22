#!/usr/bin/env python
"""Promoter deletions in the GEA's OWN results, axis by axis.

Nothing here re-estimates an association. The input is the non-SNP GEA as it stands -- the
per-record LFMM p-values in `wza_in_clq09_tile` for the 22 climate axes -- and the GEA's own
significance call, reproduced exactly from build_gea_pool.py: per class x axis, MAF > 0.05 and
a block, Bonferroni at 0.05 / n_records, a record counting as a hit on an axis if it clears
that bar in any of the three class scans (sv, smallindel, and the pooled nonsnp).

Region per record comes from results/gea_universe_regions.csv.gz: the SAME classifier the hits
carry (dissection/screen_sig_blocks.classify), run on every tested record, with its own 1 kb
promoter window and with the burden model's 2 kb one.

Per axis, three views of whether promoter deletions carry more climate signal than other
deletions, each stratified by size bin x MAF bin (Mantel-Haenszel), with 95% CIs from
resampling the GEA's own LD blocks (both hits and non-hits move with their block):

  hits    OR of being a GEA Bonferroni hit, promoter vs non-promoter deletions
  tail    the same at p < 1e-3 -- far more records per axis, so far more power than
          Bonferroni alone, still read straight off the GEA
  lambda  the ratio of genomic inflation (median chi2 / 0.455) of promoter deletions to other
          deletions: a shift of the WHOLE p distribution, not only its tail

The same three for insertions, and the hit OR for CDS insertions, as comparisons. Direction
is not tested: the saved GEA output carries p-values only, no effect sign.

Axes are correlated (axis_clusters.py), so the 22 rows are not 22 independent tests; the
per-axis Bonferroni column divides by 22 all the same.

Writes results/gea_region_by_axis_{1kb,2kb}.csv (one per promoter window). env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.abspath(os.path.join(HERE, "..", ".."))
LFMM = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/wza_in_clq09_tile"
OUT = f"{HERE}/results"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1", "pc2", "pc3"]
KEY = ["chrom", "pos", "ref_len", "alt_len"]
REPS = 1000
SIZE_BINS = [0, 1, 9, 49, 499, np.inf]
MAF_BINS = [0.05, 0.1, 0.2, 0.3, 0.5]


def axis_table(ax, regions, prom_col):
    """One row per tested non-SNP record on this axis, with the GEA's own hit call."""
    N = pd.read_csv(f"{LFMM}/lfmm_nonsnp_gen9_{ax}.csv")
    N = N[(N.MAF > 0.05) & N.block.notna() & (N.block != "")].copy()
    N["hit"] = N.pval < 0.05 / len(N)
    N["mafr"] = N.MAF.round(6)
    for cls in ("sv", "smallindel"):
        C = pd.read_csv(f"{LFMM}/lfmm_{cls}_gen9_{ax}.csv")
        C = C[(C.MAF > 0.05) & C.block.notna() & (C.block != "")]
        s = C[C.pval < 0.05 / len(C)].assign(mafr=lambda d: d.MAF.round(6))[KEY + ["mafr"]]
        s["_h"] = True
        N = N.merge(s.drop_duplicates(), on=KEY + ["mafr"], how="left")
        N["hit"] = N.hit | N._h.fillna(False).astype(bool)
        N = N.drop(columns="_h")
    N = N[N.ref_len != N.alt_len].copy()                    # MNPs: neither deletion nor insertion
    N = N.merge(regions[["chrom", "pos", "ref_len", prom_col]], on=["chrom", "pos", "ref_len"], how="left")
    N["deletion"] = N.ref_len > N.alt_len
    N["promoter"] = N[prom_col].eq("3_promoter")
    N["cds"] = N[prom_col].eq("1_CDS")
    N["tail"] = N.pval < 1e-3
    size = (N.ref_len - N.alt_len).abs()
    N["stratum"] = (pd.cut(size, SIZE_BINS, labels=False).astype(str) + "|"
                    + pd.cut(N.MAF, MAF_BINS, labels=False, include_lowest=True).astype(str))
    return N


def cell_counts(d, outcome, x):
    """block x stratum x {outcome&x, outcome&~x, ~outcome&x, ~outcome&~x} counts."""
    bl, bi = np.unique(d.block.values, return_inverse=True)
    st, si = np.unique(d.stratum.values, return_inverse=True)
    cell = (2 * (~d[outcome].values).astype(int) + (~d[x].values).astype(int))
    A = np.zeros((len(bl), len(st), 4))
    np.add.at(A, (bi, si, cell), 1)
    return A


def mh(T):
    a, b, c, e = T[..., 0], T[..., 1], T[..., 2], T[..., 3]; n = a + b + c + e
    ok = n > 0
    num = np.where(ok, a * e / np.where(ok, n, 1), 0).sum(-1)
    den = np.where(ok, b * c / np.where(ok, n, 1), 0).sum(-1)
    return num / np.where(den > 0, den, np.nan)


def or_with_ci(d, outcome, x, rng):
    A = cell_counts(d, outcome, x)
    est = mh(A.sum(0))
    # no observed outcome&x records (e.g. zero CDS-insertion hits on an axis): every bootstrap
    # draw is 0 and the two-sided p collapses to 0 -- "nothing to test", not a finding
    if A[..., 0].sum() == 0:
        return est, np.nan, np.nan, np.nan
    nb = A.shape[0]; boot = []
    for _ in range(REPS // 100):                            # batches keep the weight matrix small
        W = np.stack([np.bincount(rng.integers(0, nb, nb), minlength=nb)
                      for _ in range(100)]).astype(float)
        boot.append(mh(np.einsum("rb,bsk->rsk", W, A)))
    boot = np.concatenate(boot)
    boot = boot[np.isfinite(boot)]
    p = min(1.0, 2 * min((boot <= 1).mean(), (boot >= 1).mean())) if len(boot) else np.nan
    return est, np.percentile(boot, 2.5), np.percentile(boot, 97.5), p


def lam(p):
    return np.median(chi2.isf(np.clip(p, 1e-300, 1), 1)) / chi2.ppf(0.5, 1)


def lambda_ratio(d, rng):
    """lambda(promoter deletions) / lambda(other deletions), block bootstrap."""
    est = lam(d[d.promoter].pval.values) / lam(d[~d.promoter].pval.values)
    blocks = d.block.values; ub, bi = np.unique(blocks, return_inverse=True)
    order = np.argsort(bi); starts = np.searchsorted(bi[order], np.arange(len(ub)))
    ends = np.r_[starts[1:], len(order)]
    pv, pr = d.pval.values[order], d.promoter.values[order]
    boot = []
    for _ in range(100):                                    # medians are slow; 100 is enough for a CI
        pick = rng.integers(0, len(ub), len(ub))
        idx = np.concatenate([np.arange(starts[k], ends[k]) for k in pick])
        boot.append(lam(pv[idx][pr[idx]]) / lam(pv[idx][~pr[idx]]))
    boot = np.array(boot)
    return est, np.percentile(boot, 2.5), np.percentile(boot, 97.5)


def main(prom_col):
    regions = pd.read_csv(f"{OUT}/gea_universe_regions.csv.gz")
    rng = np.random.default_rng(1)
    rows = []
    for ax in AXES:
        N = axis_table(ax, regions, prom_col)
        Dl, In = N[N.deletion], N[~N.deletion]
        r = dict(axis=ax, n_del=len(Dl), del_hits=int(Dl.hit.sum()),
                 prom_del_hits=int((Dl.hit & Dl.promoter).sum()),
                 prom_share_hits=Dl[Dl.hit].promoter.mean() if Dl.hit.any() else np.nan,
                 prom_share_all=Dl.promoter.mean())
        for tag, d, out in [("del_hit", Dl, "hit"), ("del_tail", Dl, "tail"),
                            ("ins_hit", In, "hit"), ("ins_tail", In, "tail")]:
            if d[out].sum() < 5:
                continue
            e, lo, hi, p = or_with_ci(d, out, "promoter", rng)
            r.update({f"{tag}_or": e, f"{tag}_lo": lo, f"{tag}_hi": hi, f"{tag}_p": p})
        if In.hit.sum() >= 5:
            e, lo, hi, p = or_with_ci(In, "hit", "cds", rng)
            r.update(cdsins_hit_or=e, cdsins_hit_lo=lo, cdsins_hit_hi=hi, cdsins_hit_p=p)
        e, lo, hi = lambda_ratio(Dl, rng)
        r.update(del_lambda_ratio=e, del_lambda_lo=lo, del_lambda_hi=hi)
        rows.append(r)
        print(f"{ax:5s} del hits {r['del_hits']:>4} (promoter {r['prom_del_hits']:>3})  "
              f"hit OR {r.get('del_hit_or', np.nan):.2f} [{r.get('del_hit_lo', np.nan):.2f}-"
              f"{r.get('del_hit_hi', np.nan):.2f}]  tail OR {r.get('del_tail_or', np.nan):.2f} "
              f"[{r.get('del_tail_lo', np.nan):.2f}-{r.get('del_tail_hi', np.nan):.2f}]  "
              f"lambda ratio {e:.3f} [{lo:.3f}-{hi:.3f}]", flush=True)
    R = pd.DataFrame(rows)
    for c in [c for c in R.columns if c.endswith("_p")]:
        R[c + "_bonf22"] = (R[c] * 22).clip(upper=1)
    tag = prom_col.split("_")[1]
    R.to_csv(f"{OUT}/gea_region_by_axis_{tag}.csv", index=False)
    print(f"\nwrote {OUT}/gea_region_by_axis_{tag}.csv")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tier_1kb")
