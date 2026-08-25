#!/usr/bin/env python3
"""Quantitative hapFIRE-vs-kMate comparison on Ruth's individuals (Chr1).

hapFIRE's per-founder `_ecotype_frequency.txt` is its analog of kMate's founder
mixture `h`. Both are a probability distribution over the SAME 231 founders, so we
align by founder ID and compare with DISTRIBUTION-distance metrics (NOT Pearson
r / R² — these vectors are mostly zeros with a few large values, so correlation is
dominated by the shared tail and meaninglessly inflated):
  - effective #founders (1/Σf²) from each method
  - total variation distance  TVD = 0.5*Σ|h_kMate - h_hapFIRE|  (0=identical,
    1=disjoint; = fraction of founder mass that disagrees)
  - Hellinger distance        H = sqrt(0.5*Σ(sqrt(p)-sqrt(q))^2)  in [0,1]
  - whether the top founder agrees

hapFIRE here is Chr1-only (the hapFIRE-ready panel VCF is Chr1), so we compare to
kMate's Chr1 h (out/<sample>_Chr1.h_per_chrom.npz).

Usage: compare_ecotype_vs_h.py [--out compare_ecotype_vs_h]
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = Path("/global/scratch/users/tbellg/kmate")
KMATE_OUT = HERE.parent / "out"
KMER_PA_META = ROOT / "data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1.meta.npz"


def kmate_h(sample):
    p = KMATE_OUT / f"{sample}_Chr1.h_per_chrom.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=True)
    h = np.asarray(z["Chr1"], float)
    f = np.asarray(np.load(KMER_PA_META, allow_pickle=True)["founders"]).astype(str)
    h = h / h.sum()
    return dict(zip(f, h))


def hapfire_freq(sample):
    p = HERE / "results" / f"{sample}_ecotype_frequency.txt"
    if not p.exists():
        return None
    d = pd.read_csv(p, sep="\t", header=None, names=["founder", "freq"], dtype={0: str})
    f = d.set_index("founder")["freq"]
    f = f / f.sum()
    return f.to_dict()


def eff_n(vec):
    a = np.array(list(vec.values()), float)
    return 1.0 / np.sum(a ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "compare_ecotype_vs_h"))
    a = ap.parse_args()

    man = pd.read_csv(HERE / "samples.tsv", sep="\t")
    rows, pairs = [], {}
    for _, r in man.iterrows():
        s, lab = r["sample_id"], r["label"]
        hk, hf = kmate_h(s), hapfire_freq(s)
        if hk is None or hf is None:
            print(f"  [skip] {s}: kMate={hk is not None} hapFIRE={hf is not None}")
            continue
        founders = sorted(set(hk) | set(hf))
        k = np.array([hk.get(x, 0.0) for x in founders])
        f = np.array([hf.get(x, 0.0) for x in founders])
        k = k / k.sum(); f = f / f.sum()
        tvd = 0.5 * np.abs(k - f).sum()
        hell = np.sqrt(0.5 * np.sum((np.sqrt(k) - np.sqrt(f)) ** 2))
        topk = founders[int(np.argmax(k))]
        topf = founders[int(np.argmax(f))]
        rows.append(dict(sample=s, label=lab,
                         kMate_effN=eff_n(hk), hapFIRE_effN=eff_n(hf),
                         TVD=tvd, hellinger=hell, top_founder_agree=(topk == topf),
                         kMate_top=topk, hapFIRE_top=topf))
        pairs[s] = (k, f, lab)

    if not rows:
        print("No samples with BOTH outputs yet.")
        return
    df = pd.DataFrame(rows).sort_values(["label", "sample"]).reset_index(drop=True)
    pd.set_option("display.width", 160)
    print(df.to_string(index=False))
    df.to_csv(a.out + ".csv", index=False)
    print(f"\nsaved -> {a.out}.csv")

    # scatter h vs ecotype_freq, one panel per sample
    n = len(pairs)
    cols = min(n, 4); rows_ = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows_, cols, figsize=(3.4 * cols, 3.4 * rows_), squeeze=False)
    cmap = {"inbred": "#4C72B0", "outcrossed": "#DD8452"}
    for ax, (s, (k, f, lab)) in zip(axes.flat, pairs.items()):
        ax.scatter(f, k, s=8, alpha=0.5, color=cmap[lab], edgecolor="none")
        m = max(k.max(), f.max()) * 1.05
        ax.plot([0, m], [0, m], ls="--", c="0.5", lw=1)
        tvd = 0.5 * np.abs(k - f).sum()
        ax.set_title(f"{s.split('_')[0].replace('MEAJM003-','')} ({lab[:3]})  TVD={tvd:.3f}",
                     fontsize=9)
        ax.set_xlabel("hapFIRE ecotype freq"); ax.set_ylabel("kMate h")
    for ax in axes.flat[len(pairs):]:
        ax.axis("off")
    fig.suptitle("Per-founder frequency: kMate h vs hapFIRE ecotype_frequency (Chr1)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(a.out + "_scatter.png", dpi=150, bbox_inches="tight")
    print(f"saved -> {a.out}_scatter.png")

    # eff_n agreement
    fig2, ax = plt.subplots(figsize=(7, 4))
    x = np.arange(len(df)); w = 0.38
    ax.bar(x - w/2, df.kMate_effN, w, label="kMate (h)", color="#55A868")
    ax.bar(x + w/2, df.hapFIRE_effN, w, label="hapFIRE (ecotype)", color="#C44E52")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s.split('_')[0].replace('MEAJM003-','')}\n{l[:3]}"
                        for s, l in zip(df["sample"], df.label)], fontsize=8)
    ax.set_ylabel("effective # founders (1/Σf²)")
    ax.set_title("Effective #founders: kMate vs hapFIRE")
    ax.legend(fontsize=9)
    fig2.tight_layout(); fig2.savefig(a.out + "_effn.png", dpi=150, bbox_inches="tight")
    print(f"saved -> {a.out}_effn.png")


if __name__ == "__main__":
    main()
