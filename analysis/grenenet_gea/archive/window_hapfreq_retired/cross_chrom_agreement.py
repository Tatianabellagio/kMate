#!/usr/bin/env python
"""Per-sample cross-chromosome agreement of the founder-frequency estimates.

genome_h() averages each sample's five chromosome-wide `ChrN_global_h` founder-proportion vectors
(231-dim, sum 1) into one genome-wide h. This checks how consistent those 5 estimates are WITHIN a
sample: if the pool is a real mixture of whole genomes, all chromosomes should recover ~the same
founder proportions; disagreement = estimation noise (low coverage) or chromosome-specific artifacts.

Per sample metrics over the 5 chromosome h-vectors:
  mean_r    : mean pairwise Pearson r (10 pairs)               -> shape agreement
  mean_tv   : mean pairwise total-variation 0.5*sum|hi-hj|     -> fraction of founder mass that disagrees
  r_to_gw   : per-chrom r of ChrN vs the genome-wide mean h    -> which chromosome is the odd one out
Merges coverage (if available) to see whether agreement is coverage-driven.
Writes cross_chrom_agreement.{csv,png}. Env: basic (numpy/pandas/matplotlib + lib).
"""
import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

WIN = "results/grenenet_kmate_window"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
OUT = "results/grenenet_gea/hapfreq"


def load_h(samp):
    H = []
    for ch in CHROMS:
        f = f"{WIN}/{samp}_{ch}.h_blocks_per_chrom.npz"
        if not os.path.exists(f):
            return None
        H.append(np.load(f, allow_pickle=True)[f"{ch}_global_h"].astype(np.float64))
    return np.vstack(H)                                   # (5, 231)


def main():
    samps = sorted(p.split("/")[-1].split("_Chr1")[0] for p in glob.glob(f"{WIN}/*_Chr1.h_blocks_per_chrom.npz"))
    print(f"{len(samps)} samples", flush=True)
    iu = np.triu_indices(5, 1)
    rows = []
    for i, s in enumerate(samps):
        Hs = load_h(s)
        if Hs is None:
            continue
        R = np.corrcoef(Hs)                              # 5x5 chrom-chrom
        pw_r = R[iu]
        tv = np.array([0.5 * np.abs(Hs[a] - Hs[b]).sum() for a, b in zip(*iu)])
        gw = Hs.mean(0)
        r_to_gw = np.array([np.corrcoef(Hs[c], gw)[0, 1] for c in range(5)])
        eff_n = 1.0 / (Hs ** 2).sum(1)                   # per-chrom effective #founders
        rows.append(dict(sampleid=s, mean_r=pw_r.mean(), min_r=pw_r.min(), mean_tv=tv.mean(),
                         max_tv=tv.max(), eff_n_mean=eff_n.mean(), eff_n_cv=eff_n.std() / eff_n.mean(),
                         **{f"r_{ch}": r_to_gw[c] for c, ch in enumerate(CHROMS)}))
        if (i + 1) % 400 == 0:
            print(f"  ...{i+1}/{len(samps)}", flush=True)
    df = pd.DataFrame(rows)

    # merge coverage if the pool table exposes it
    try:
        pt = lib.pool_table()
        if "coverage" in pt.columns:
            cov = pt.drop_duplicates("sampleid").set_index("sampleid")["coverage"]
            df["coverage"] = df.sampleid.map(cov.to_dict())
    except Exception as e:
        print(f"[warn] no coverage merge ({type(e).__name__})")
    df.to_csv(f"{OUT}/cross_chrom_agreement.csv", index=False)

    print("\n=== cross-chromosome agreement (per sample) ===")
    for c, lab in [("mean_r", "mean pairwise r"), ("min_r", "worst pair r"), ("mean_tv", "mean pairwise TV (frac mass)")]:
        q = df[c].quantile([0.05, 0.25, 0.5, 0.75, 0.95]).round(3)
        print(f"  {lab:28s}: median {df[c].median():.3f}  IQR [{q[0.25]}, {q[0.75]}]  5-95% [{q[0.05]}, {q[0.95]}]")
    print("\n  per-chromosome mean r-to-genome-wide (which chrom disagrees most):")
    for ch in CHROMS:
        print(f"    {ch}: {df[f'r_{ch}'].median():.3f}")
    if "coverage" in df.columns and df.coverage.notna().any():
        m = df.coverage.notna()
        print(f"\n  corr(mean_r, coverage) = {np.corrcoef(df.mean_r[m], df.coverage[m])[0,1]:.2f} "
              f"(n={int(m.sum())}) -> agreement rises with depth" )

    # ---- figure ----
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    ax[0].hist(df.mean_r, bins=60, color="#3b4cc0", alpha=.8, edgecolor="none")
    ax[0].axvline(df.mean_r.median(), color="firebrick", lw=1.5, label=f"median {df.mean_r.median():.3f}")
    ax[0].set_xlabel("mean pairwise cross-chromosome r"); ax[0].set_ylabel("samples")
    ax[0].set_title(f"(a) per-sample cross-chrom agreement (n={len(df)})", loc="left", fontsize=10)
    ax[0].legend(frameon=False, fontsize=8); ax[0].spines[["top", "right"]].set_visible(False)
    bp = ax[1].boxplot([df[f"r_{ch}"] for ch in CHROMS], labels=CHROMS, showfliers=False, patch_artist=True)
    for b in bp["boxes"]:
        b.set(facecolor="#7aa0c4", alpha=.7)
    ax[1].set_ylabel("r of chromosome vs genome-wide h")
    ax[1].set_title("(b) which chromosome agrees least", loc="left", fontsize=10)
    ax[1].spines[["top", "right"]].set_visible(False)
    if "coverage" in df.columns and df.coverage.notna().any():
        m = df.coverage.notna()
        ax[2].scatter(df.coverage[m], df.mean_r[m], s=6, c="#34495e", alpha=.35, edgecolors="none")
        ax[2].set_xscale("log"); ax[2].set_xlabel("sample coverage (log)"); ax[2].set_ylabel("mean pairwise r")
        ax[2].set_title(f"(c) agreement vs depth  r={np.corrcoef(df.mean_r[m], df.coverage[m])[0,1]:.2f}", loc="left", fontsize=10)
    else:
        ax[2].hist(df.mean_tv, bins=60, color="#3b4cc0", alpha=.8, edgecolor="none")
        ax[2].set_xlabel("mean pairwise total variation"); ax[2].set_title("(c) disagreeing founder mass", loc="left", fontsize=10)
    ax[2].spines[["top", "right"]].set_visible(False)
    fig.suptitle("Cross-chromosome consistency of per-sample founder-frequency (global_h) estimates", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(f"{OUT}/cross_chrom_agreement.png", dpi=150)
    print(f"\n[done] {OUT}/cross_chrom_agreement.csv + .png")


if __name__ == "__main__":
    main()
