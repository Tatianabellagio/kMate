#!/usr/bin/env python
"""Promoter deletions and garden temperature, asked at the SV EVENT level.

Why a second test. The GEA -- and deletion_promoter_enrichment.py, which inherits its unit --
tests every arch3 RECORD separately. arch3 splits one SV into one record per assembly path,
so a single insertion can appear as many low-frequency "alleles": the signal is split,
fragments fall below MAF 0.05 and are never tested, and insertions fragment far more than
deletions, so a record-level deletion-vs-insertion comparison is tilted before it starts.
The purging work met exactly this: 70% of record-level "cold-specific" insertions were
fragments of shared events (notebooks/sv_purging_fixed.ipynb). The burden model that found
the promoter-deletion x temperature effect (notebooks/sv_climate_maf.ipynb: +0.133% per
variant per degC, permutation p < 0.001; one-stage LRT p = 0.0046) uses EVENTS -- the
alleles of a pangenome bubble collapsed with `truvari collapse` defaults.

So this uses the same unit and the same data as that work, read-only from
r1_sv_negative_selection/results/sv_adaptive/:
  garden_sgslope_event.npz   per variant x garden selection coefficient s (garden-level
                             logit-frequency slope across generations; "the one to use"),
                             col 0 = p0; events for SVs, a 10% record sample of short indels
  sv_event_catalog.csv       event position, size, carriers, carrier-origin climate, context
Filters are theirs (called in >= 90% of founders, >= 2 carriers). Context codes 0..5 =
CDS, UTR, intron, TE, promoter, intergenic (first match wins, as in the burden model).

Per variant: beta = OLS slope of s on garden bio1 across the 31 gardens (> 0 = the ALT allele
-- for a deletion event, the deletion -- does relatively better in warmer gardens).

  A  THE BURDEN CONTRAST, per variant. Within deletions: is beta higher for promoter events
     than intergenic ones? Same contrast as the burden model (promoter vs intergenic within
     a class), so it is free of reference polarisation -- both are deletions relative to
     Col-0. Adjusted for carrier-origin climate (a variant carried by warm-origin founders
     rises in warm gardens by hitchhiking on its lineage), logit p0, log carriers and log
     size. CIs by resampling 1 Mb genomic windows (nearby events share haplotypes).
     The same for insertions, and for all six contexts against intergenic.
  B  ENRICHMENT AMONG THE MOST CLIMATE-SELECTED EVENTS. "Hits" = events whose slope t-statistic
     is in the top 5% by |t|; is promoter over-represented among deletion hits? MH over
     carrier-climate x p0 x size strata, window-bootstrap CI, and the same for insertions.

Writes results/event_promoter_deletion_{contrast,enrichment}.csv. env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd
import statsmodels.formula.api as smf

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = f"{HERE}/results"
SVA = os.path.abspath(os.path.join(HERE, "..", "..", "r1_sv_negative_selection", "results", "sv_adaptive"))
CTX = ["CDS", "UTR", "intron", "TE", "promoter", "intergenic"]
WINDOW = 1_000_000
REPS = 2000


def slopes(S, x):
    """Row-wise OLS slope and t of S (variants x gardens) on x, NaN-aware."""
    ok = np.isfinite(S)
    n = ok.sum(1).astype(float)
    X = np.where(ok, x[None, :], np.nan)
    xm = np.nanmean(X, 1, keepdims=True); ym = np.nanmean(np.where(ok, S, np.nan), 1, keepdims=True)
    dx = np.where(ok, X - xm, 0.0); dy = np.where(ok, S - ym, 0.0)
    sxx = (dx ** 2).sum(1); b = (dx * dy).sum(1) / np.where(sxx > 0, sxx, np.nan)
    res = np.where(ok, dy - b[:, None] * dx, 0.0)
    se = np.sqrt((res ** 2).sum(1) / np.maximum(n - 2, 1) / np.where(sxx > 0, sxx, np.nan))
    b[n < 20] = np.nan
    return b, b / se


def load():
    z = np.load(f"{SVA}/garden_sgslope_event.npz", allow_pickle=True)
    b, t = slopes(z["p"][:, 1:].astype(float), z["bio1"].astype(float))
    D = pd.DataFrame(dict(kind=z["kind"].astype(str), ctx=[CTX[c] for c in z["ctx"]],
                          n=z["n"].astype(int), b1=z["b1"].astype(float),
                          origin=z["origin"].astype(str), p0=z["p"][:, 0].astype(float),
                          cat_row=z["cat_row"].astype(int), beta=b, t=t))
    C = pd.read_csv(f"{SVA}/sv_event_catalog.csv", usecols=["chrom", "pos", "size_bp"])
    ev = D.cat_row >= 0
    D.loc[ev, "chrom"] = C.chrom.values[D.cat_row[ev]]
    D.loc[ev, "pos"] = C.pos.values[D.cat_row[ev]]
    D.loc[ev, "size_bp"] = C.size_bp.values[D.cat_row[ev]]
    D = D[ev & np.isfinite(D.beta)].copy()
    D["window"] = D.chrom + ":" + (D.pos // WINDOW).astype(int).astype(str)
    p = D.p0.clip(1e-3, 1 - 1e-3)
    D["lp0"] = np.log(p / (1 - p)); D["lcar"] = np.log10(D.n); D["lsize"] = np.log10(D.size_bp)
    return D


def contrast(D, kind, shared_only):
    """beta ~ context (vs intergenic) + carrier climate + logit p0 + log carriers + log size."""
    d = D[(D.kind == kind) & (D.origin.eq("shared") if shared_only else True)].copy()
    d["ctx"] = pd.Categorical(d.ctx, categories=["intergenic"] + [c for c in CTX if c != "intergenic"])
    f = "beta ~ C(ctx) + b1 + lp0 + lcar + lsize"
    fit = smf.ols(f, data=d).fit()
    terms = [f"C(ctx)[T.{c}]" for c in CTX if c != "intergenic"]
    wins = d.window.unique(); groups = {w: np.where(d.window.values == w)[0] for w in wins}
    rng = np.random.default_rng(1); boot = []
    for _ in range(REPS):
        idx = np.concatenate([groups[w] for w in rng.choice(wins, len(wins), replace=True)])
        try:
            boot.append(smf.ols(f, data=d.iloc[idx]).fit().params.reindex(terms).values)
        except Exception:                                          # noqa: BLE001
            continue
    boot = np.array(boot)
    rows = []
    for j, tm in enumerate(terms):
        bb = boot[:, j]; bb = bb[np.isfinite(bb)]
        rows.append(dict(kind=kind, shared_only=shared_only, context=tm.split(".")[1][:-1],
                         n=int((d.ctx == tm.split(".")[1][:-1]).sum()),
                         n_intergenic=int((d.ctx == "intergenic").sum()),
                         coef=fit.params[tm], ci_lo=np.percentile(bb, 2.5),
                         ci_hi=np.percentile(bb, 97.5),
                         boot_p=min(1.0, 2 * min((bb <= 0).mean(), (bb >= 0).mean())),
                         b1_coef=fit.params["b1"]))
    return rows


def enrichment(D, kind, top=0.05):
    d = D[D.kind == kind].copy()
    thr = np.nanquantile(np.abs(D[D.kind.isin(["del", "ins"])].t), 1 - top)
    d["hit"] = np.abs(d.t) >= thr
    d["x"] = d.ctx.eq("promoter")
    d["stratum"] = (pd.qcut(d.b1, 3, labels=False, duplicates="drop").astype(str) + "|"
                    + pd.qcut(d.p0, 4, labels=False, duplicates="drop").astype(str) + "|"
                    + pd.cut(d.size_bp, [0, 100, 500, 2000, np.inf], labels=False).astype(str))

    def mh(dd, w=None):
        w = np.ones(len(dd)) if w is None else w
        num = den = 0.0
        for s, g in dd.assign(_w=w).groupby("stratum"):
            a = g._w[g.hit & g.x].sum(); b = g._w[g.hit & ~g.x].sum()
            c = g._w[~g.hit & g.x].sum(); e = g._w[~g.hit & ~g.x].sum(); n = a + b + c + e
            if n:
                num += a * e / n; den += b * c / n
        return num / den if den else np.nan

    orr = mh(d)
    wins = d.window.unique(); pos = {w: i for i, w in enumerate(wins)}
    wi = d.window.map(pos).values; rng = np.random.default_rng(2); boot = []
    for _ in range(REPS):
        cnt = np.bincount(rng.integers(0, len(wins), len(wins)), minlength=len(wins))
        boot.append(mh(d, cnt[wi].astype(float)))
    boot = np.array(boot); boot = boot[np.isfinite(boot)]
    h = d[d.hit]
    return dict(kind=kind, top=top, n=len(d), n_hits=int(d.hit.sum()),
                promoter_frac_hits=h.x.mean(), promoter_frac_nonhits=d[~d.hit].x.mean(),
                mh_or=orr, ci_lo=np.percentile(boot, 2.5), ci_hi=np.percentile(boot, 97.5),
                boot_p=min(1.0, 2 * min((boot <= 1).mean(), (boot >= 1).mean())),
                promoter_hits_warm_favoured=(h[h.x].beta > 0).mean() if h.x.any() else np.nan,
                other_hits_warm_favoured=(h[~h.x].beta > 0).mean())


def main():
    D = load()
    print(f"{len(D):,} SV events with a slope: {dict(D.kind.value_counts())}; "
          f"{D.window.nunique()} 1-Mb windows", flush=True)
    rows = []
    for kind in ("del", "ins"):
        for sh in (False, True):
            rows += contrast(D, kind, sh)
    A = pd.DataFrame(rows)
    B = pd.DataFrame([enrichment(D, k, top) for k in ("del", "ins") for top in (0.05, 0.01)])
    A.to_csv(f"{OUT}/event_promoter_deletion_contrast.csv", index=False)
    B.to_csv(f"{OUT}/event_promoter_deletion_enrichment.csv", index=False)
    pd.set_option("display.width", 240)
    print("\nA  beta (slope of s on garden bio1) by context vs intergenic, adjusted for carrier climate,"
          " logit p0, carriers, size; 1-Mb window bootstrap")
    print(A.round(4).to_string(index=False))
    print("\nB  promoter share among the most climate-selected events (top |t|), MH over "
          "carrier-climate x p0 x size strata")
    print(B.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
