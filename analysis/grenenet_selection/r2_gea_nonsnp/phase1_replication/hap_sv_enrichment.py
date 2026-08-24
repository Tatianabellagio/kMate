#!/usr/bin/env python
"""Step 3: are climate-associated clq0.9 haploblocks enriched for SVs?

Joins the Kendall-tau haploblock GEA (kendall_hap_gen9_bio1.csv) to the per-block SV
tags (hap_gen9.registry.csv), then asks whether haploblocks that associate with climate
are over-represented for SV-containing blocks.

Confounds handled (all matter here because global-h projection inflates significance
genome-wide, and SV blocks differ in geometry from SNP-only blocks):
  * BLOCK LEVEL, not haplotype: a big multi-haplotype block is ONE observation.
  * WITHIN-BLOCK MULTIPLE TESTING: SV blocks carry more haplotypes, so their min-p is
    lower just from more draws. Block p = Sidak(min_p, n_hap) removes that bias before
    anything else. (This is the dominant confound; do not skip it.)
  * COVARIATE CONTROL: SV blocks are longer / higher n_eff / carry more haplotypes;
    the has_sv effect is reported net of MAF + log(block_len) + n_eff + log(n_hap).

Three readouts (agreement = a real effect):
  (1) THRESHOLD-FREE  linreg  -log10(p_block) ~ has_sv + covariates  (+ Mann-Whitney)
  (2) THRESHOLD Fisher OR   sig(BH q<0.05) x has_sv
  (3) THRESHOLD logistic    sig ~ has_sv + covariates  (adjusted has_sv OR)
Env: kmate (statsmodels optional; the linreg/logistic degrade gracefully).

  python hap_sv_enrichment.py            # after run_kendall.py --class hap --gen 9
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import lib

CM = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/class_matrices"
KEN = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/kendall/kendall_hap_gen9_bio1.csv"
OUT = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/hap_sv_enrichment.csv"


def bh(p):
    p = np.asarray(p, float); m = len(p); o = p.argsort(); q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def fisher(sig, has_sv):
    a = int((sig & has_sv).sum()); b = int((sig & ~has_sv).sum())
    c = int((~sig & has_sv).sum()); d = int((~sig & ~has_sv).sum())
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    return dict(sv_sig=a, snp_sig=b, sv_ns=c, snp_ns=d, OR=float(orr), p=float(p),
                sig_rate_sv=a / max(a + c, 1), sig_rate_snp=b / max(b + d, 1))


def size_matched(blk, seed=0):
    """SIZE-MATCHED null (the correct test; a length-covariate logistic OVER-CONTROLS because
    has_sv and block length are near-collinear -- SV blocks live in a longer size regime, so
    regression extrapolates across non-overlapping support and fakes a depletion).

    Match each SV block to its nearest unused SNP-only block by (log10 len, log10 nvar); compare
    %BH-sig and median -log10(p_block); paired Wilcoxon + McNemar. Also a within-length-decile
    permutation p (shuffle has_sv within deciles) as a second, coarser check."""
    rng = np.random.default_rng(seed)
    b = blk.copy()
    b["logL"] = np.log10((b.block_end - b.block_start).clip(lower=1))
    b["logV"] = np.log10(b.unit_nvar.clip(lower=1))
    b["nlp"] = -np.log10(b.p_block.clip(1e-300, 1)); b["sig"] = (b.q_block < 0.05).astype(int)
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
    nd = int((sv_sig != mm_sig).sum()); mc = stats.binomtest(int((sv_sig > mm_sig).sum()), nd) if nd else None
    # within-length-decile permutation on mean nlp difference
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


def main():
    ap = __import__("argparse").ArgumentParser()
    ap.add_argument("--gea", default=KEN, help="per-haploblock GEA csv (needs hap_id,pval); "
                    "default=Kendall. Pass the LFMM output for the structure-corrected version.")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    globals()["OUT"] = a.out
    ken = pd.read_csv(a.gea)
    reg = pd.read_csv(f"{CM}/hap_gen9.registry.csv")
    d = ken.merge(reg, on="hap_id", suffixes=("", "_reg"))
    d = d[np.isfinite(d.pval)].copy()
    print(f"GEA input: {a.gea}")
    lam = np.median(stats.chi2.isf(d.pval.clip(1e-300, 1), 1)) / stats.chi2.ppf(.5, 1)
    print(f"haploblock tests: {len(d):,} haps over {d.block.nunique():,} blocks | "
          f"lambda_GC={lam:.2f} | raw p<0.05: {(d.pval<0.05).mean():.1%}")

    # ---- BLOCK level: min-p haplotype per block + within-block Sidak(n_hap) ----
    g = d.groupby("block")
    lead = d.loc[g.pval.idxmin()].copy()                 # most climate-associated hap per block
    lead["n_hap"] = g.size().reindex(lead.block).to_numpy()
    lead["p_block"] = 1 - np.power(1 - lead.pval.clip(0, 1), lead.n_hap)   # Sidak within block
    lead["q_block"] = bh(lead.p_block.to_numpy())
    blk = lead.reset_index(drop=True)
    print(f"blocks: {len(blk):,} | has_sv {int(blk.has_sv.sum()):,} ({blk.has_sv.mean():.1%})")

    # (1) UNCONTROLLED Fisher OR -- reported ONLY to show the raw hint (confounded by block size)
    f = fisher((blk.q_block < 0.05).to_numpy(), blk.has_sv.to_numpy())
    print(f"\n(1) RAW (uncontrolled, confounded by block size): SV %BHsig {f['sig_rate_sv']:.1%} vs "
          f"SNP {f['sig_rate_snp']:.1%} | Fisher OR={f['OR']:.2f} p={f['p']:.2g}")
    print("    -> SV blocks are 9.6x longer, so this raw OR is a geometry artifact, NOT an SV effect.")

    # (2) SIZE-MATCHED null -- the correct test (a length-covariate logistic OVER-controls;
    #     see size_matched docstring). This is the number to trust.
    sm = size_matched(blk)
    print(f"\n(2) SIZE-MATCHED null ({sm['n_pairs']:,} pairs, median |dlogLen|={sm['match_dlogL']:.3f}):")
    print(f"    %BHsig     SV={sm['sig_rate_sv']:.1%}  matched-SNP={sm['sig_rate_match']:.1%}")
    print(f"    med -log10p SV={sm['med_nlp_sv']:.3f}  matched-SNP={sm['med_nlp_match']:.3f} "
          f"| paired Wilcoxon p={sm['wilcoxon_p']:.2g} | McNemar p={sm['mcnemar_p']:.2g}")
    print(f"    within-length-decile permutation: mean Delta={sm['perm_meanDiff']:+.3f} p={sm['perm_p']:.3f}")
    verdict = ("marginal/unstable hint" if (sm['wilcoxon_p'] < 0.05) != (sm['perm_p'] < 0.05)
               else "enriched" if sm['perm_p'] < 0.05 and sm['med_nlp_sv'] > sm['med_nlp_match']
               else "no enrichment")
    print(f"    VERDICT: {verdict} of SV climate-association at matched block size "
          f"(on a lambda={lam:.1f}-inflated scan).")

    rows = [dict(lambda_gc=lam, n_blocks=len(blk), raw_fisher_OR=f['OR'], raw_fisher_p=f['p'],
                 raw_sig_rate_sv=f['sig_rate_sv'], raw_sig_rate_snp=f['sig_rate_snp'],
                 verdict=verdict, **{f"matched_{k}": v for k, v in sm.items()})]
    pd.DataFrame(rows).to_csv(OUT, index=False)
    blk.to_csv(OUT.replace(".csv", "_blocks.csv"), index=False)
    print(f"\n[done] {OUT} (+ _blocks.csv)")
    top = blk[blk.has_sv].sort_values("p_block").head(15)
    cols = ["block", "chrom", "block_start", "block_end", "maf", "unit_n_eff", "n_hap",
            "n_sv", "sv_bp_max"] + (["tau"] if "tau" in blk.columns else []) + ["pval", "p_block", "q_block"]
    print("\nTop climate-associated SV-containing haploblocks:")
    print(top[cols].to_string(index=False))


if __name__ == "__main__":
    main()
