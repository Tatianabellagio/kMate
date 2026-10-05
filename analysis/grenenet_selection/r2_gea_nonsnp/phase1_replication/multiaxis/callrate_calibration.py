#!/usr/bin/env python
"""Is the GEA calibrated the same way for well-called and poorly-called records?

The panel is 78 founders genotyped from assemblies + 153 from short reads by PanGenie, so a
record's allele frequency is `AF = (h @ var_pa) / (h @ var_called)` over however many founders
are called there. Records called in fewer founders carry less information AND -- because
records sharing a called-set share the denominator exactly -- share a common error direction
across pools. 6.1% of records share one single pattern (the 78 assembly founders).

That is a latent covariate, and the risk is that such records are systematically more likely
to clear a threshold built from all records at once. FLC is a worked instance: filling its 153
missing founders from the phase-1 panel drops its garden-level climate gradient from
r = -0.668 to -0.574 and makes kMate agree with hapFIRE at r = 0.997.

This measures the size of the effect genome-wide, the only way that settles whether the
single per-(class x axis) Bonferroni is safe: lambda, tail inflation and hit rate computed
WITHIN call-rate strata, and within MAF strata as well, since call rate and frequency are
correlated and MAF is the known confounder.

Strata: fully called (>=0.99) | 0.90-0.99 | 0.50-0.90 | the 78-only set (0.3377 +- 0.002) |
other < 0.50.

Out: results/multiaxis/callrate_calibration.csv  (one row per class x axis x call stratum x
MAF stratum), plus ..._bymaf.csv marginalised over axes.

env: kmate. Compute node.
"""
from __future__ import annotations
import argparse
import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..")))
import lib                                                        # noqa: E402

MULTI = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis"
WZAIN = f"{MULTI}/_superseded_pre_cr50/wza_in_clq09_tile"   # the PRE-floor run, on purpose: this script measures what the 0.50 call-rate floor removes
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1", "pc2", "pc3"]
MAF_BINS = [0.05, 0.10, 0.20, 0.50001]
QS = [0.5, 0.9, 0.99, 0.999]
SUB78 = 78 / 231


def call_rates(kind: str) -> pd.DataFrame:
    """(chrom,pos,ref_len,alt_len) -> founders called / 231, for every af_store record."""
    idx = np.load(f"{lib.AF_STORE}/index_{kind}.npz", allow_pickle=True)
    ch = idx["chrom"].astype(str)
    rl = idx["ref_len"].astype(np.int64)
    al = idx["alt_len"].astype(np.int64)
    n_cal = np.zeros(len(ch), np.int64)
    off = 0
    for c in ("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"):
        n = int((ch == c).sum())
        base = f"{lib.PROJ}/panel/arch3/{c.lower()}/var_pa_231_arch3_{c.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        prl = meta["ref_len"].astype(np.int64)
        pal = meta["alt_len"].astype(np.int64)
        issnp = (prl == 1) & (pal == 1)
        panel = np.where(issnp if kind == "snp" else ~issnp)[0]
        assert panel.size == n, f"{c}/{kind}: af_store {n} vs panel {panel.size}"
        assert np.array_equal(prl[panel], rl[off:off + n]), f"{c}/{kind}: record order differs"
        vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        n_cal[off:off + n] = np.asarray((vc[:, panel] > 0).sum(0)).ravel()
        off += n
        print(f"  call rates {kind} {c}", flush=True)
    return pd.DataFrame(dict(chrom=ch, pos=idx["pos"], ref_len=rl, alt_len=al,
                             call_rate=n_cal / 231.0))


def stratum(cr: np.ndarray) -> np.ndarray:
    s = np.full(len(cr), "other <0.50", dtype=object)
    s[cr < 0.50] = "other <0.50"
    s[np.isclose(cr, SUB78, atol=0.002)] = "the 78-only set"
    s[(cr >= 0.50) & (cr < 0.90)] = "0.50-0.90"
    s[(cr >= 0.90) & (cr < 0.99)] = "0.90-0.99"
    s[cr >= 0.99] = "fully called"
    return s


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", default="nonsnp,sv,snp")
    a = ap.parse_args()
    rows = []
    for cls in a.classes.split(","):
        kind = "snp" if cls == "snp" else "nonsnp"
        CR = call_rates(kind)
        for ax in AXES:
            f = f"{WZAIN}/lfmm_{cls}_gen9_{ax}.csv"
            if not os.path.exists(f):
                continue
            d = pd.read_csv(f, usecols=["chrom", "pos", "ref_len", "alt_len", "MAF", "pval"])
            d = d[d.MAF > 0.05]
            thr = 0.05 / len(d)                       # the production threshold, all records
            d = d.merge(CR, on=["chrom", "pos", "ref_len", "alt_len"], how="left")
            d = d[d.call_rate.notna()]
            d["strat"] = stratum(d.call_rate.values)
            d["mb"] = pd.cut(d.MAF, MAF_BINS, right=False)
            d["c2"] = chi2.isf(np.clip(d.pval.values, 1e-300, 1.0), 1)
            for (st, mb), g in d.groupby(["strat", "mb"], observed=True):
                if len(g) < 2000:
                    continue
                r = dict(cls=cls, axis=ax, strat=st, maf=f"{mb.left:.2f}-{mb.right:.2f}",
                         n=len(g), med_call_rate=float(g.call_rate.median()),
                         lam=float(np.median(g.c2) / 0.4549364),
                         hits=int((g.pval < thr).sum()),
                         hit_per_10k=1e4 * float((g.pval < thr).mean()))
                for q in QS:
                    r[f"q{q}"] = float(np.quantile(g.c2, q) / chi2.ppf(q, 1))
                rows.append(r)
            print(f"  {cls:7s} {ax:5s} done", flush=True)
    D = pd.DataFrame(rows)
    D.to_csv(f"{MULTI}/callrate_calibration.csv", index=False)

    pd.set_option("display.width", 240)
    ORDER = ["fully called", "0.90-0.99", "0.50-0.90", "the 78-only set", "other <0.50"]
    for cls in D.cls.unique():
        S = D[D.cls == cls]
        w = S.groupby(["maf", "strat"]).apply(
            lambda x: pd.Series({"n": x.n.sum(),
                                 "lam": np.average(x.lam, weights=x.n),
                                 "q0.999": np.average(x["q0.999"], weights=x.n),
                                 "hit_per_10k": 1e4 * x.hits.sum() / x.n.sum()}),
            include_groups=False).reset_index()
        print(f"\n===== {cls}: calibration within MAF x call-rate strata "
              f"(lambda, tail at the 99.9th pct, Bonferroni hit rate)\n")
        for m in sorted(w.maf.unique()):
            t = w[w.maf == m].set_index("strat").reindex(ORDER).dropna(how="all")
            print(f"  MAF {m}"); print(t.round(3).to_string(), "\n")
    print(f"wrote {MULTI}/callrate_calibration.csv")


if __name__ == "__main__":
    main()
