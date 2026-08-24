#!/usr/bin/env python
"""Cross-site synthesis of the per-site temporal SV-enrichment (clq0.9 blocks).

Reads every analysis/grenenet_selection/r1_sv_negative_selection/results/site_temporal/site{N}_sv_enrichment.json produced by
site_sv_enrichment.py and answers the headline question across sites:
  Is there SV enrichment in the temporally-selected haploblocks -- anywhere, consistently?

Two readouts, per site + pooled:
  B (region enrichment): size-matched fold + empirical p of SV-bearing blocks among the
    top-1% selected clq0.9 blocks. Pooled sign-test / count of sites with fold>1 & p<0.05.
  A (per-variant effect): p0-matched median|s| SV/SNP ratio (is it consistently >1?).

Outputs: site_temporal/cross_site_enrichment.csv, .png (figure), and a printed verdict.
Run in the `basic` env (matplotlib works there; plotting env's is broken).
  PY=<basic env>; $PY aggregate_sites_enrichment.py
"""
from __future__ import annotations
import glob, json, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

ST = f"{lib.GEA}/r1_sv_negative_selection/results/site_temporal"


def main():
    rows = []
    for f in sorted(glob.glob(f"{ST}/site*_sv_enrichment.json")):
        d = json.load(open(f))
        rows.append(dict(
            site=d["site"], n_selected=d["B_n_selected"], n_scored=d["B_n_scored"],
            hasSV_sel=d["B_hasSV_selected"], hasSV_null=d["B_hasSV_null_median"],
            hasSV_nonsel=d["B_hasSV_nonselected"], fold=d["B_fold"], p_emp=d["B_p_emp"],
            raw_OR=d["B_raw_OR_hasSV"], raw_p=d["B_raw_fisher_p"],
            A_ratio=d["A_p0matched_ratio_sv_snp"],
            s_snp=d["A_median_abs_s"]["snp"], s_sv=d["A_median_abs_s"]["sv"]))
    df = pd.DataFrame(rows)
    clim = lib.load_climate()["bio1"]
    df["bio1"] = df.site.map(clim).astype(float)
    # trajectory length per site (from pool metas): 3=gen1-3, 2, 1=Δp-only
    ng = {}
    for g in (1, 2, 3):
        m = pd.read_csv(f"{lib.GEA}/common/results/pool_matrices/pool_gen{g}_nonsnp.meta.csv")
        for s in m.site.unique():
            ng[int(s)] = ng.get(int(s), 0) + 1
    df["n_gens"] = df.site.map(ng)
    df = df.sort_values("bio1").reset_index(drop=True)

    pd.set_option("display.width", 160)
    print(df[["site", "bio1", "n_gens", "n_selected", "hasSV_sel", "hasSV_null",
              "fold", "p_emp", "raw_OR", "A_ratio"]].to_string(index=False,
              float_format=lambda x: f"{x:.3f}"))

    n = len(df)
    n_fold_gt1 = int((df.fold > 1).sum())
    n_sig = int((df.p_emp < 0.05).sum())
    n_sig_up = int(((df.p_emp < 0.05) & (df.fold > 1)).sum())
    from scipy.stats import binomtest
    sign_p = binomtest(n_fold_gt1, n, 0.5).pvalue
    print(f"\n=== VERDICT (region enrichment, {n} sites) ===")
    print(f"sites with fold>1 (SV-enriched direction): {n_fold_gt1}/{n} "
          f"(sign-test vs 0.5: p={sign_p:.3f})")
    print(f"sites with p_emp<0.05: {n_sig}  |  enriched & significant (fold>1 & p<0.05): {n_sig_up}")
    print(f"fold: median={df.fold.median():.3f}  mean={df.fold.mean():.3f}  "
          f"range=[{df.fold.min():.2f}, {df.fold.max():.2f}]")
    print(f"=== Part A (per-variant |s| SV/SNP, p0-matched) ===")
    print(f"ratio>1 in {int((df.A_ratio>1).sum())}/{n} sites  |  "
          f"median={df.A_ratio.median():.3f}  range=[{df.A_ratio.min():.2f}, {df.A_ratio.max():.2f}]")
    df.to_csv(f"{ST}/cross_site_enrichment.csv", index=False)

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    sig = df.p_emp < 0.05
    ax[0].axhline(1, color="grey", ls="--", lw=1)
    ax[0].scatter(df.bio1, df.fold, c=np.where(sig, "#ee6677", "#4477aa"), s=40)
    for _, r in df.iterrows():
        ax[0].annotate(int(r.site), (r.bio1, r.fold), fontsize=7,
                       xytext=(2, 2), textcoords="offset points")
    ax[0].set_xlabel("site bio1 (mean annual T, °C)"); ax[0].set_ylabel("SV-bearing fold (size-matched)")
    ax[0].set_title("B. region SV-enrichment per site\n(red = p_emp<0.05)")
    ax[1].axhline(1, color="grey", ls="--", lw=1)
    ax[1].scatter(df.bio1, df.A_ratio, c="#228833", s=40)
    ax[1].set_xlabel("site bio1 (°C)"); ax[1].set_ylabel("per-variant |s| SV/SNP (p0-matched)")
    ax[1].set_title("A. SVs move more than SNPs?")
    ax[2].hist(df.fold, bins=np.linspace(0.4, 1.6, 25), color="#bbbbbb")
    ax[2].axvline(1, color="k", ls="--", lw=1)
    ax[2].axvline(df.fold.median(), color="#ee6677", lw=2, label=f"median {df.fold.median():.2f}")
    ax[2].set_xlabel("SV-bearing fold across sites"); ax[2].set_ylabel("n sites")
    ax[2].set_title("C. fold distribution"); ax[2].legend()
    fig.tight_layout(); fig.savefig(f"{ST}/cross_site_enrichment.png", dpi=130)
    print(f"\n-> {ST}/cross_site_enrichment.csv + .png")


if __name__ == "__main__":
    main()
