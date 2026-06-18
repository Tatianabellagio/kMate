#!/usr/bin/env python
"""Full cap-vs-no-cap comparison across ALL blocks (not just CAM5).

Does the SNP cap actually change the genome-wide WZA result, or only avoid the NaN
windows? Compares 4 regimes on gen1 SNP x bio1 (Kendall):
  deg2_nocap (+floor), deg2_cap2000, deg7_nocap (phase-1), deg7_cap2000
For each: NaN count, GIF (genomic inflation), #Bonferroni, #BH-sig(q<.05), and
between-regime agreement: Spearman of -log10 p, BH-sig set overlap (Jaccard), and
status-flippers (BH-sig in one but not the other).
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd
from scipy.stats import norm, spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wza_core as wc

KEN = ("/global/scratch/users/tbellg/kmate/results/grenenet_gea/"
       "phase1_replication/kendall/kendall_snp_gen1_bio1.csv")
OUT = "/global/scratch/users/tbellg/kmate/results/grenenet_gea/wza_investigation"


def bh(p):
    """Benjamini-Hochberg q-values for a p-value Series (NaN-safe)."""
    q = pd.Series(np.nan, index=p.index)
    v = p.dropna().sort_values()
    n = len(v)
    raw = v.to_numpy() * n / np.arange(1, n + 1)
    q.loc[v.index] = np.minimum.accumulate(raw[::-1])[::-1]
    return q


def gif(p):
    """Genomic inflation factor for one-sided p (z^2 ~ chi2_1; median/0.4549)."""
    z = norm.ppf(1 - p.dropna().clip(1e-300, 1 - 1e-16))
    return np.median(z**2) / 0.4549549


def main():
    df = pd.read_csv(KEN)
    df = df[df.MAF >= 0.05].copy()
    base = wc.raw_wza(df, cap=None)
    cap2k = wc.raw_wza(df, cap=2000)

    regimes = {
        "deg2_nocap_floor": (base, 2, True),
        "deg2_cap2000":     (cap2k, 2, False),
        "deg7_nocap":       (base, 7, False),
        "deg7_cap2000":     (cap2k, 7, False),
    }
    res = {}
    summ = []
    for name, (rw, deg, floor) in regimes.items():
        w, _ = wc.apply_correction(rw, deg=deg, sd_floor=floor)
        w["q"] = bh(w["Z_pVal"])
        res[name] = w.set_index("block")
        p = w["Z_pVal"]
        summ.append(dict(regime=name, blocks=len(w), nan=int(p.isnull().sum()),
                         GIF=round(gif(p), 3),
                         bonferroni=int((p < 0.05/len(p.dropna())).sum()),
                         BH_q05=int((w["q"] < 0.05).sum())))
    summ = pd.DataFrame(summ)
    summ.to_csv(f"{OUT}/cap_compare_summary.csv", index=False)
    print("=== PER-REGIME SUMMARY ===")
    print(summ.to_string(index=False))

    # --- agreement between regimes (focus cap vs nocap at each degree) ---
    print("\n=== AGREEMENT: cap vs no-cap (same degree) ===")
    pairs = [("deg2_nocap_floor", "deg2_cap2000"),
             ("deg7_nocap", "deg7_cap2000"),
             ("deg2_cap2000", "deg7_cap2000")]
    agree = []
    for a, b in pairs:
        A, B = res[a], res[b]
        idx = A.index.intersection(B.index)
        pa, pb = A.loc[idx, "Z_pVal"], B.loc[idx, "Z_pVal"]
        m = pa.notna() & pb.notna()
        rho = spearmanr(-np.log10(pa[m].clip(1e-300)), -np.log10(pb[m].clip(1e-300))).statistic
        sa = set(A.index[A["q"] < 0.05]); sb = set(B.index[B["q"] < 0.05])
        jac = len(sa & sb) / max(len(sa | sb), 1)
        agree.append(dict(regimeA=a, regimeB=b, spearman_logp=round(rho, 4),
                          BHsig_A=len(sa), BHsig_B=len(sb), shared=len(sa & sb),
                          only_A=len(sa - sb), only_B=len(sb - sa), jaccard=round(jac, 3)))
    agree = pd.DataFrame(agree)
    agree.to_csv(f"{OUT}/cap_compare_agreement.csv", index=False)
    print(agree.to_string(index=False))

    # --- the NaN windows that no-cap (no floor) silently drops ---
    w_nofloor, _ = wc.apply_correction(base, deg=2, sd_floor=False)
    nanb = w_nofloor[w_nofloor["Z_pVal"].isnull()]
    print(f"\n=== {len(nanb)} blocks get NO p-value under deg2 no-cap (no floor) ===")
    print(nanb[["block", "SNPs_raw", "Z"]].sort_values("SNPs_raw", ascending=False).head(20).to_string(index=False))
    nanb.to_csv(f"{OUT}/cap_compare_nan_blocks.csv", index=False)

    # --- top-50 candidate list stability (deg2): which blocks enter/leave with cap ---
    A, B = res["deg2_nocap_floor"], res["deg2_cap2000"]
    topA = set(A["Z_pVal"].nsmallest(50).index); topB = set(B["Z_pVal"].nsmallest(50).index)
    print(f"\n=== deg2 top-50 blocks: cap vs nocap — shared {len(topA & topB)}/50, "
          f"only-nocap {len(topA - topB)}, only-cap {len(topB - topA)} ===")

    _plot(res)
    print(f"\nwrote outputs to {OUT}")


def _plot(res):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 6))
    for a, (deg, A, B) in zip(ax, [("deg2", "deg2_nocap_floor", "deg2_cap2000"),
                                    ("deg7", "deg7_nocap", "deg7_cap2000")]):
        pass
    specs = [("deg-2", "deg2_nocap_floor", "deg2_cap2000"),
             ("deg-7", "deg7_nocap", "deg7_cap2000")]
    for a, (lab, an, bn) in zip(ax, specs):
        A, B = res[an], res[bn]
        idx = A.index.intersection(B.index)
        x = -np.log10(A.loc[idx, "Z_pVal"].clip(1e-300))
        y = -np.log10(B.loc[idx, "Z_pVal"].clip(1e-300))
        snp = A.loc[idx, "SNPs_raw"]
        sc = a.scatter(x, y, c=np.log10(snp), s=8, alpha=.5, cmap="viridis")
        lim = max(x.max(), y.max())
        a.plot([0, lim], [0, lim], "r--", lw=1)
        a.set(xlabel=f"-log10 p  NO CAP ({lab})", ylabel=f"-log10 p  CAP2000 ({lab})",
              title=f"{lab}: per-block p, cap vs no-cap")
        plt.colorbar(sc, ax=a, label="log10(SNPs/block)")
    fig.suptitle("Full cap-vs-no-cap comparison (gen1 SNP, bio1) — color = block SNP count")
    fig.tight_layout(); fig.savefig(f"{OUT}/fig8_cap_vs_nocap_allblocks.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
