#!/usr/bin/env python
"""Are climate-associated clq0.9 blocks (kendall/lfmm/binomial -> WZA, STATUS_clq90.md)
enriched for SVs?

Tests at the SAME unit the association was tested on: the clq0.9 block IS the WZA
window, so (unlike hap_sv_enrichment.py, which had multiple haplotypes per block and
needed a within-block Sidak correction) each block already carries exactly one p-value
per (model, class=nonsnp) WZA run. No representative-choice step needed.

Composition (has_sv / n_sv / n_indel / block span) comes from
driver_passenger/block_composition_kendall_gen9_bio1.csv, which tags ALL 58,308 clq0.9
blocks regardless of which WZA run scored them (it's built from the raw per-class
records, not from any model's p-values) -- confirmed to cover 100% of the 42,677
nonsnp-class blocks.

SV blocks are far larger (more variants -> lower min-p just from more draws, same
confound as hap_sv_enrichment.py), so:
  (1) RAW Fisher OR sig(BH q<0.05) x has_sv        -- shown only to demonstrate the
      confound, NOT trusted.
  (2) SIZE-MATCHED null (the number to trust): nearest-neighbour match each SV block to
      an unused non-SV block by (log10 block length, log10 n_var); compare %BH-sig and
      median -log10p; paired Wilcoxon + McNemar; within-length-decile permutation.

Runs across all 3 models (kendall/lfmm/binomial) x deg2 (canonical, primary).

Env: kmate.
  python sv_enrichment_clq90.py                 # deg2, all 3 models
  python sv_enrichment_clq90.py --deg deg7cap2000
"""
from __future__ import annotations
import argparse, glob, os, re, sys
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import lib

WZA_DIR = f"{lib.GEA}/phase1_replication/results/clq90/wza"
COMPOSITION = f"{lib.GEA}/driver_passenger/results/block_composition_kendall_gen9_bio1.csv"
OUT = f"{lib.GEA}/phase1_replication/results/clq90/sv_enrichment_clq90.csv"


def bh(p):
    p = np.asarray(p, float); m = len(p); o = p.argsort(); q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def fisher(sig, has_sv):
    a = int((sig & has_sv).sum()); b = int((sig & ~has_sv).sum())
    c = int((~sig & has_sv).sum()); d = int((~sig & ~has_sv).sum())
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    return dict(sv_sig=a, nonsv_sig=b, sv_ns=c, nonsv_ns=d, OR=float(orr), p=float(p),
                sig_rate_sv=a / max(a + c, 1), sig_rate_nonsv=b / max(b + d, 1))


def size_matched(blk, seed=0):
    """Nearest-neighbour size-matched null (mirrors hap_sv_enrichment.py::size_matched).
    A length-covariate logistic over-controls because has_sv is collinear with block
    length/n_var (SV blocks live in a longer/denser regime) -- match instead."""
    rng = np.random.default_rng(seed)
    b = blk.copy()
    b["nlp"] = -np.log10(b.pval.clip(1e-300, 1)); b["sig"] = (b.q < 0.05).astype(int)
    sv = b.has_sv.to_numpy()
    A = b[sv]; B = b[~sv].reset_index(drop=True); Bxy = np.column_stack([B.logL, B.logV])
    used = np.zeros(len(B), bool); ai, bj = [], []
    for idx in A.index.to_numpy()[A.nlp.to_numpy().argsort()]:
        dd = ((Bxy - np.array([b.logL[idx], b.logV[idx]])) ** 2).sum(1); dd[used] = np.inf
        j = int(dd.argmin()); used[j] = True; ai.append(idx); bj.append(B.index[j])
    sv_nlp = b.loc[ai, "nlp"].to_numpy(); mm_nlp = B.loc[bj, "nlp"].to_numpy()
    sv_sig = b.loc[ai, "sig"].to_numpy(); mm_sig = B.loc[bj, "sig"].to_numpy()
    dL = float(np.median(np.abs(b.loc[ai, "logL"].to_numpy() - B.loc[bj, "logL"].to_numpy())))
    w = stats.wilcoxon(sv_nlp, mm_nlp, zero_method="zsplit")
    nd = int((sv_sig != mm_sig).sum())
    mc = stats.binomtest(int((sv_sig > mm_sig).sum()), nd) if nd else None
    strata = pd.qcut(b.logL, 10, labels=False, duplicates="drop").to_numpy()

    def stat(mask):
        vals = [b.nlp[mask & (strata == s)].mean() - b.nlp[~mask & (strata == s)].mean()
                for s in np.unique(strata) if (mask & (strata == s)).any() and (~mask & (strata == s)).any()]
        return float(np.mean(vals))

    obs = stat(sv); null = np.empty(2000)
    for k in range(2000):
        perm = np.zeros(len(b), bool)
        for s in np.unique(strata):
            si = np.where(strata == s)[0]; perm[rng.choice(si, int(sv[si].sum()), replace=False)] = True
        null[k] = stat(perm)
    return dict(n_pairs=len(ai), match_dlogL=dL,
                sig_rate_sv=float(sv_sig.mean()), sig_rate_match=float(mm_sig.mean()),
                med_nlp_sv=float(np.median(sv_nlp)), med_nlp_match=float(np.median(mm_nlp)),
                wilcoxon_p=float(w.pvalue), mcnemar_p=float(mc.pvalue) if mc else float("nan"),
                perm_meanDiff=obs, perm_p=float((np.abs(null) >= abs(obs)).mean()))


def run_one(model: str, deg: str, comp: pd.DataFrame) -> dict:
    path = f"{WZA_DIR}/wza_{model}_nonsnp_gen9_bio1_{deg}.csv"
    w = pd.read_csv(path).rename(columns={"index": "block", "Z_pVal": "pval", "SNPs": "n_var"})
    w = w[np.isfinite(w.pval)].copy()
    w["q"] = bh(w.pval.to_numpy())
    d = w.merge(comp[["block", "has_sv", "n_sv", "n_indel", "start", "end"]], on="block", how="inner")
    assert len(d) == len(w), f"{model}/{deg}: {len(w) - len(d)} blocks failed to match composition"
    d["logL"] = np.log10((d.end - d.start).clip(lower=1))
    d["logV"] = np.log10(d.n_var.clip(lower=1))

    f = fisher((d.q < 0.05).to_numpy(), d.has_sv.to_numpy())
    sm = size_matched(d)
    verdict = ("marginal/unstable hint" if (sm['wilcoxon_p'] < 0.05) != (sm['perm_p'] < 0.05)
               else "enriched" if sm['perm_p'] < 0.05 and sm['med_nlp_sv'] > sm['med_nlp_match']
               else "no enrichment")
    print(f"\n=== model={model} deg={deg} ===  n_blocks={len(d):,} has_sv={int(d.has_sv.sum()):,} "
          f"({d.has_sv.mean():.1%})  BH-sig={int((d.q < 0.05).sum())}")
    print(f"  (1) RAW: sig_rate SV={f['sig_rate_sv']:.1%} nonSV={f['sig_rate_nonsv']:.1%} "
          f"OR={f['OR']:.2f} p={f['p']:.2g}  [confounded by block size -- not trusted]")
    print(f"  (2) SIZE-MATCHED ({sm['n_pairs']:,} pairs, median |dlogLen|={sm['match_dlogL']:.3f}): "
          f"sig SV={sm['sig_rate_sv']:.1%} matched={sm['sig_rate_match']:.1%} | "
          f"med -log10p SV={sm['med_nlp_sv']:.3f} matched={sm['med_nlp_match']:.3f} | "
          f"Wilcoxon p={sm['wilcoxon_p']:.2g} McNemar p={sm['mcnemar_p']:.2g} perm p={sm['perm_p']:.3f}")
    print(f"  VERDICT: {verdict}")

    return dict(model=model, deg=deg, n_blocks=len(d), n_has_sv=int(d.has_sv.sum()),
                n_bh_sig=int((d.q < 0.05).sum()), raw_OR=f['OR'], raw_p=f['p'],
                raw_sig_rate_sv=f['sig_rate_sv'], raw_sig_rate_nonsv=f['sig_rate_nonsv'],
                verdict=verdict, **{f"matched_{k}": v for k, v in sm.items()})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", nargs="+", default=["kendall", "lfmm", "binomial"])
    ap.add_argument("--deg", default="deg2", choices=["deg2", "deg7cap2000", "both"])
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    degs = ["deg2", "deg7cap2000"] if a.deg == "both" else [a.deg]

    comp = pd.read_csv(COMPOSITION)
    rows = [run_one(m, d, comp) for d in degs for m in a.models]
    out = pd.DataFrame(rows)
    out.to_csv(a.out, index=False)
    print(f"\n[done] {a.out}")


if __name__ == "__main__":
    main()
