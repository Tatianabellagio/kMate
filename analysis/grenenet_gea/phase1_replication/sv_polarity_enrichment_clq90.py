#!/usr/bin/env python
"""Follow-up to sv_enrichment_clq90.py: split has_sv by insertion/deletion polarity.

Motivation: sv_enrichment_clq90.py found a discordant signal for lfmm/binomial
(permutation test on the *mean* -log10p shift flags something; McNemar, the direct
test of excess BH-sig crossing, does not) -- exactly the small-tail/bulk-shift pattern
`SV_TEMPORAL_PURGING_SUMMARY.md` found per-variant, where the real effect lives in
**SV insertions only** (climate-graded purifying selection), deletions null.

This asks the same question at the BLOCK level: is the has_sv WZA-enrichment hint
actually an insertion-only effect? Blocks are split into three DISJOINT SV categories
(ins-only / del-only / both), each tested against the same non-SV control pool
(has_sv==False) with the identical size-matched null as before.

Env: kmate.
  python sv_polarity_enrichment_clq90.py
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import lib

WZA_DIR = f"{lib.GEA}/phase1_replication/results/clq90/wza"
COMPOSITION = f"{lib.GEA}/driver_passenger/results/block_composition_kendall_gen9_bio1.csv"
POLARITY = f"{lib.GEA}/phase1_replication/results/clq90/block_polarity_composition.csv"
OUT = f"{lib.GEA}/phase1_replication/results/clq90/sv_polarity_enrichment_clq90.csv"


def bh(p):
    p = np.asarray(p, float); m = len(p); o = p.argsort(); q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def size_matched(b: pd.DataFrame, treat: np.ndarray, ctrl: np.ndarray, seed=0):
    """Nearest-neighbour size-matched null: match each treat block to the nearest
    unused ctrl block by (logL, logV). Mirrors sv_enrichment_clq90.py::size_matched
    but with an explicit control pool (not just ~treat)."""
    rng = np.random.default_rng(seed)
    A = b[treat]; B = b[ctrl].reset_index(drop=True)
    if len(A) < 5 or len(B) < 5:
        return dict(n_pairs=len(A), match_dlogL=float("nan"), sig_rate_sv=float("nan"),
                    sig_rate_match=float("nan"), med_nlp_sv=float("nan"), med_nlp_match=float("nan"),
                    wilcoxon_p=float("nan"), mcnemar_p=float("nan"), perm_meanDiff=float("nan"),
                    perm_p=float("nan"))
    Bxy = np.column_stack([B.logL, B.logV]); used = np.zeros(len(B), bool); ai, bj = [], []
    for idx in A.index.to_numpy()[A.nlp.to_numpy().argsort()]:
        dd = ((Bxy - np.array([b.logL[idx], b.logV[idx]])) ** 2).sum(1); dd[used] = np.inf
        j = int(dd.argmin()); used[j] = True; ai.append(idx); bj.append(B.index[j])
    sv_nlp = b.loc[ai, "nlp"].to_numpy(); mm_nlp = B.loc[bj, "nlp"].to_numpy()
    sv_sig = b.loc[ai, "sig"].to_numpy(); mm_sig = B.loc[bj, "sig"].to_numpy()
    dL = float(np.median(np.abs(b.loc[ai, "logL"].to_numpy() - B.loc[bj, "logL"].to_numpy())))
    w = stats.wilcoxon(sv_nlp, mm_nlp, zero_method="zsplit")
    nd = int((sv_sig != mm_sig).sum())
    mc = stats.binomtest(int((sv_sig > mm_sig).sum()), nd) if nd else None
    # permutation null: stratify treat-vs-ctrl pool by logL decile, shuffle treat labels within strata
    sub = b[treat | ctrl].copy()
    strata = pd.qcut(sub.logL, 10, labels=False, duplicates="drop").to_numpy()

    def stat(mask):
        vals = [sub.nlp[mask & (strata == s)].mean() - sub.nlp[~mask & (strata == s)].mean()
                for s in np.unique(strata) if (mask & (strata == s)).any() and (~mask & (strata == s)).any()]
        return float(np.mean(vals)) if vals else float("nan")

    tsub = np.asarray(sub.index.isin(A.index))
    obs = stat(tsub); null = np.empty(2000)
    for k in range(2000):
        perm = np.zeros(len(sub), bool)
        for s in np.unique(strata):
            si = np.where(strata == s)[0]; perm[rng.choice(si, int(tsub[si].sum()), replace=False)] = True
        null[k] = stat(perm)
    perm_p = float((np.abs(null) >= abs(obs)).mean()) if np.isfinite(obs) else float("nan")
    return dict(n_pairs=len(ai), match_dlogL=dL,
                sig_rate_sv=float(sv_sig.mean()), sig_rate_match=float(mm_sig.mean()),
                med_nlp_sv=float(np.median(sv_nlp)), med_nlp_match=float(np.median(mm_nlp)),
                wilcoxon_p=float(w.pvalue), mcnemar_p=float(mc.pvalue) if mc else float("nan"),
                perm_meanDiff=obs, perm_p=perm_p)


def run_one(model: str, deg: str, comp: pd.DataFrame, pol: pd.DataFrame, cls: str = "nonsnp") -> list[dict]:
    path = f"{WZA_DIR}/wza_{model}_{cls}_gen9_bio1_{deg}.csv"
    w = pd.read_csv(path).rename(columns={"index": "block", "Z_pVal": "pval", "SNPs": "n_var"})
    w = w[np.isfinite(w.pval)].copy()
    w["q"] = bh(w.pval.to_numpy())
    d = w.merge(comp[["block", "has_sv", "start", "end"]], on="block", how="inner")
    d = d.merge(pol, on="block", how="left")
    d[["n_ins_sv", "n_del_sv"]] = d[["n_ins_sv", "n_del_sv"]].fillna(0)
    assert len(d) == len(w)
    d["logL"] = np.log10((d.end - d.start).clip(lower=1))
    d["logV"] = np.log10(d.n_var.clip(lower=1))
    d["nlp"] = -np.log10(d.pval.clip(1e-300, 1))
    d["sig"] = (d.q < 0.05).astype(int)

    ins_only = (d.n_ins_sv > 0) & (d.n_del_sv == 0)
    del_only = (d.n_del_sv > 0) & (d.n_ins_sv == 0)
    both = (d.n_ins_sv > 0) & (d.n_del_sv > 0)
    control = ~d.has_sv  # no SV at all (may still contain indels)

    rows = []
    for label, mask in [("ins_only_sv", ins_only.to_numpy()),
                         ("del_only_sv", del_only.to_numpy()),
                         ("both_sv", both.to_numpy())]:
        sm = size_matched(d, mask, control.to_numpy())
        n_bh = int(d.sig[mask].sum())
        print(f"[{model}/{deg}/{cls}] {label}: n_blocks={int(mask.sum()):,} BH-sig={n_bh} | "
              f"size-matched sig SV={sm['sig_rate_sv']:.1%} matched={sm['sig_rate_match']:.1%} "
              f"McNemar p={sm['mcnemar_p']:.2g} Wilcoxon p={sm['wilcoxon_p']:.2g} perm p={sm['perm_p']:.3f}")
        rows.append(dict(model=model, deg=deg, cls=cls, category=label, n_blocks=int(mask.sum()),
                          n_bh_sig=n_bh, **{f"matched_{k}": v for k, v in sm.items()}))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="+", default=["kendall", "lfmm", "binomial"])
    ap.add_argument("--cls", nargs="+", default=["nonsnp"], choices=["nonsnp", "snp"])
    ap.add_argument("--deg", default="deg2", choices=["deg2", "deg7cap2000", "both"])
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    degs = ["deg2", "deg7cap2000"] if a.deg == "both" else [a.deg]

    comp = pd.read_csv(COMPOSITION)
    pol = pd.read_csv(POLARITY)
    rows = []
    for d in degs:
        for m in a.models:
            for c in a.cls:
                print(f"\n=== model={m} deg={d} cls={c} ===")
                rows += run_one(m, d, comp, pol, cls=c)
    out = pd.DataFrame(rows)
    out.to_csv(a.out, index=False)
    print(f"\n[done] {a.out}")


if __name__ == "__main__":
    main()
