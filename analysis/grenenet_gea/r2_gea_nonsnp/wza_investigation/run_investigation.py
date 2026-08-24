#!/usr/bin/env python
"""WZA investigation on GrENE-Net phase-1 — gen1 SNPs, bio1.

Answers three worries:
  1. (doc, not here) LD threshold of the 16,674 blocks = hapFIRE BigLD, independent
     r2=0.1 + fine CLQcut r2=0.5 (density). See README.md.
  2. (doc) RepAdapt windows = genes +500bp flank, with a per-gene SNP cap. See README.md.
  3. Does WZA over-weight low-SNP windows (e.g. CAM5, ~16 SNPs)? Tested below by
     binning corrected p-values by SNP count + scatter of Z_pVal vs SNP count.

Also sweeps polynomial degree {2,7} x SNP cap {none, 2000, q95, q75} and reports
NaN rate, Bonferroni hits, and the CAM5 block (2_1265 gene 3'; 2_1264 pericentro null).
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except ModuleNotFoundError:
    HAVE_MPL = False

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wza_core as wc

KENDALL = ("/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/"
           "phase1_replication/kendall/kendall_snp_gen1_bio1.csv")
OUT = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/wza_investigation/results"
os.makedirs(OUT, exist_ok=True)
CAM5 = ["2_1265", "2_1264"]   # gene 3' (real hit) ; gene-midpoint pericentro null


def bonf(w):
    p = w["Z_pVal"].dropna()
    return int((p < 0.05 / len(p)).sum()) if len(p) else 0


def cam5_rows(w):
    w = w.copy()
    w["rank"] = w["Z_pVal"].rank()
    return w[w["block"].isin(CAM5)][["block", "SNPs_raw", "Z", "sd_pred",
                                     "Z_pVal", "rank"]]


def main():
    df = pd.read_csv(KENDALL)
    df = df[df["MAF"] >= 0.05].copy()   # phase-1 maf05; matches canonical WZA
    print(f"gen1 SNP records (MAF>=0.05): {len(df):,} | blocks: {df['block'].nunique():,}")

    # --- raw weighted-Z per block at several caps (computed once each) ---
    base = wc.raw_wza(df, cap=None)
    snp_q95 = int(base["SNPs_raw"].quantile(0.95))
    snp_q75 = int(base["SNPs_raw"].quantile(0.75))
    print(f"SNPs/block: median {base.SNPs_raw.median():.0f}  q75 {snp_q75}  "
          f"q95 {snp_q95}  max {base.SNPs_raw.max()}")
    caps = {"nocap": None, "cap2000": 2000, f"capq95_{snp_q95}": snp_q95,
            f"capq75_{snp_q75}": snp_q75}
    raw = {k: (base if v is None else wc.raw_wza(df, cap=v)) for k, v in caps.items()}

    # --- sweep degree x cap ---
    rows = []
    cam5_all = []
    keep_for_plots = {}
    for capname, rw in raw.items():
        for deg in (2, 7):
            w, diag = wc.apply_correction(rw, deg=deg)
            tag = f"deg{deg}_{capname}"
            nan = int(w["Z_pVal"].isnull().sum())
            negsd = int(w["neg_sd"].sum())
            rows.append(dict(variant=tag, deg=deg, cap=capname, blocks=len(w),
                             nan_pvals=nan, neg_sd_windows=negsd, bonferroni=bonf(w)))
            c = cam5_rows(w); c.insert(0, "variant", tag); cam5_all.append(c)
            w.to_csv(f"{OUT}/wza_{tag}.csv", index=False)
            if capname == "nocap" or (capname == "cap2000" and deg == 2):
                keep_for_plots[tag] = (w, diag)
    summary = pd.DataFrame(rows)
    summary.to_csv(f"{OUT}/sweep_summary.csv", index=False)
    cam5_df = pd.concat(cam5_all, ignore_index=True)
    cam5_df.to_csv(f"{OUT}/cam5_across_variants.csv", index=False)
    print("\n=== SWEEP SUMMARY ===\n", summary.to_string(index=False))
    print("\n=== CAM5 ACROSS VARIANTS ===\n", cam5_df.to_string(index=False))

    # --- worry #3: significance vs SNP count (overweighting test) ---
    bins = [0, 5, 10, 20, 50, 100, 200, 500, 1000, 10**9]
    lbl = ["2-5", "6-10", "11-20", "21-50", "51-100", "101-200",
           "201-500", "501-1000", ">1000"]
    ov_rows = []
    for tag in ["deg2_nocap", "deg7_nocap", "deg2_cap2000"]:
        w = pd.read_csv(f"{OUT}/wza_{tag}.csv")
        w["bin"] = pd.cut(w["SNPs_raw"], bins=bins, labels=lbl)
        for b, sub in w.groupby("bin", observed=True):
            p = sub["Z_pVal"].dropna()
            ov_rows.append(dict(variant=tag, snp_bin=b, n=len(sub),
                                frac_p_lt_05=(p < 0.05).mean() if len(p) else np.nan,
                                frac_p_lt_001=(p < 0.001).mean() if len(p) else np.nan,
                                median_p=p.median() if len(p) else np.nan))
    ov = pd.DataFrame(ov_rows)
    ov.to_csv(f"{OUT}/overweighting_by_snpbin.csv", index=False)
    print("\n=== SIGNIFICANCE BY SNP-COUNT BIN (uniform => no SNP-count bias) ===")
    print(ov.to_string(index=False))

    if HAVE_MPL:
        _plots(base, keep_for_plots, ov, lbl)
        print(f"\nwrote outputs + figures to {OUT}")
    else:
        print(f"\n[matplotlib missing — CSVs written to {OUT}, figures skipped]")


def _plots(base, keep, ov, lbl):
    # Fig 1: SNP/block dist + SD-correction curves (where SD<=0)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
    ax[0].hist(np.log10(base["SNPs_raw"]), bins=60, color="steelblue")
    ax[0].set(xlabel="log10(SNPs per block)", ylabel="blocks",
              title=f"SNP/block (n={len(base)}, max={base.SNPs_raw.max()})")
    for tag, (w, diag) in keep.items():
        if "nocap" not in tag:
            continue
        deg = 7 if "deg7" in tag else 2
        xs = np.linspace(diag["x"].min(), base["SNPs_raw"].max(), 400)
        sd_poly = np.poly1d(np.polyfit(diag["x"], diag["sd"], deg))
        ax[1].plot(xs, sd_poly(xs), label=f"deg{deg} fit")
    ax[1].scatter(keep["deg2_nocap"][1]["x"], keep["deg2_nocap"][1]["sd"],
                  s=6, color="grey", alpha=.4, label="rolling SD (data)")
    ax[1].axhline(0, color="red", ls="--", lw=1)
    ax[1].set(xlabel="SNPs per block", ylabel="predicted SD of Z",
              title="SNP-number correction: SD curve\n(red=0 -> NaN below)")
    ax[1].legend(fontsize=8)
    # zoom tail
    for tag, (w, diag) in keep.items():
        if "nocap" not in tag:
            continue
        deg = 7 if "deg7" in tag else 2
        xs = np.linspace(diag["x"].min(), base["SNPs_raw"].max(), 400)
        sd_poly = np.poly1d(np.polyfit(diag["x"], diag["sd"], deg))
        ax[2].plot(xs, sd_poly(xs), label=f"deg{deg} fit")
    ax[2].axhline(0, color="red", ls="--", lw=1)
    ax[2].set(xlim=(200, base["SNPs_raw"].max()), ylim=(-2, 3),
              xlabel="SNPs per block (tail)", ylabel="predicted SD of Z",
              title="tail zoom: deg curves -> negative SD")
    ax[2].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig1_snp_dist_and_sd_curves.png", dpi=130)
    plt.close(fig)

    # Fig 2: Z_pVal vs SNP count (overweighting scatter)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.5), sharey=True)
    for a, tag in zip(ax, ["deg2_nocap", "deg7_nocap", "deg2_cap2000"]):
        w = pd.read_csv(f"{OUT}/wza_{tag}.csv")
        a.scatter(w["SNPs_raw"], -np.log10(w["Z_pVal"].clip(lower=1e-300)),
                  s=5, alpha=.3)
        a.set(xscale="log", xlabel="SNPs per block (log)",
              title=tag, ylabel="-log10(Z_pVal)")
        a.axhline(-np.log10(0.05 / len(w)), color="red", ls="--", lw=1)
    fig.suptitle("Per-block significance vs SNP count (worry #3)")
    fig.tight_layout(); fig.savefig(f"{OUT}/fig2_sig_vs_snpcount.png", dpi=130)
    plt.close(fig)

    # Fig 3: fraction significant per SNP bin
    fig, ax = plt.subplots(figsize=(9, 5))
    for tag in ["deg2_nocap", "deg7_nocap", "deg2_cap2000"]:
        s = ov[ov.variant == tag].set_index("snp_bin").reindex(lbl)
        ax.plot(lbl, s["frac_p_lt_05"], marker="o", label=tag)
    ax.axhline(0.05, color="grey", ls=":", label="0.05 (uniform null)")
    ax.set(xlabel="SNPs per block (bin)", ylabel="fraction Z_pVal < 0.05",
           title="Overweighting test: fraction significant by window SNP count")
    ax.legend(); plt.xticks(rotation=45)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig3_overweighting_by_snpbin.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
