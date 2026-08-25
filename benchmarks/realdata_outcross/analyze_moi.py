#!/usr/bin/env python3
"""MOI / outcrossing readout from kMate's founder-mixture h.

For each individual we load the per-chrom h vectors (global EM, saved by
per_sample_per_chrom.py as <sample>_<chrom>.h_per_chrom.npz), average them to a
genome-wide founder composition, and derive MOI statistics:

  eff_n_founders = 1 / sum(h^2)   (inverse-Simpson effective # founders)
  max_h          = largest single-founder fraction
  top2_sum       = fraction explained by the top 2 founders

Inbred  -> h ~ one-hot  -> eff_n ~ 1, max_h ~ 1.
Outcross-> h spread     -> eff_n >= 2, max_h ~ 0.5 (F1) or less.

We compare kMate's call to the hapFIRE labels in manifest.tsv and print a
per-sample table + the threshold that separates the two classes; a figure shows
eff_n_founders by label and each sample's top-founder composition.

Usage: analyze_moi.py [--dir out] [--manifest manifest.tsv] [--out moi_summary]
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]


def load_founders(kmer_pa_meta):
    m = np.load(kmer_pa_meta, allow_pickle=True)
    return np.asarray(m["founders"]).astype(str)


def genomewide_h(sample, out_dir):
    """Mean of the per-chrom h vectors (each on the founder simplex)."""
    hs = []
    for chrom in CHROMS:
        p = out_dir / f"{sample}_{chrom}.h_per_chrom.npz"
        if not p.exists():
            continue
        z = np.load(p, allow_pickle=True)
        # single-chrom run -> the chrom key holds h; fall back to first non-meta key
        key = chrom if chrom in z.files else next(k for k in z.files if k != "founders")
        hs.append(np.asarray(z[key], dtype=float))
    if not hs:
        return None, 0
    H = np.vstack(hs)
    return H.mean(axis=0), len(hs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(HERE / "out"))
    ap.add_argument("--manifest", default=str(HERE / "manifest.tsv"))
    ap.add_argument("--kmer-pa-meta",
                    default="data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1.meta.npz")
    ap.add_argument("--out", default=str(HERE / "moi_summary"))
    a = ap.parse_args()

    out_dir = Path(a.dir)
    founders = load_founders(a.kmer_pa_meta)
    man = pd.read_csv(a.manifest, sep="\t")

    rows = []
    comps = {}  # sample -> genome-wide h
    for _, r in man.iterrows():
        s = r["sample_id"]
        h, nchr = genomewide_h(s, out_dir)
        if h is None:
            print(f"  [skip] {s}: no h files yet")
            continue
        h = h / h.sum()
        order = np.argsort(h)[::-1]
        eff_n = 1.0 / np.sum(h ** 2)
        top = [(founders[i], float(h[i])) for i in order[:3]]
        rows.append(dict(
            sample=s, label=r["label"], n_chrom=nchr,
            eff_n_founders=eff_n, max_h=float(h[order[0]]),
            top2_sum=float(h[order[0]] + h[order[1]]),
            top1=f"{top[0][0]}:{top[0][1]:.2f}",
            top2=f"{top[1][0]}:{top[1][1]:.2f}",
            top3=f"{top[2][0]}:{top[2][1]:.2f}",
        ))
        comps[s] = (h, order)

    if not rows:
        print("No completed samples yet.")
        return
    df = pd.DataFrame(rows).sort_values(["label", "eff_n_founders"])
    pd.set_option("display.width", 160)
    print(df.to_string(index=False))
    csv = a.out + ".csv"
    df.to_csv(csv, index=False)
    print(f"\nsaved -> {csv}")

    # separation: best eff_n threshold between the two labels
    inb = df[df.label == "inbred"].eff_n_founders
    out = df[df.label == "outcrossed"].eff_n_founders
    if len(inb) and len(out):
        thr = (inb.max() + out.min()) / 2
        sep = inb.max() < out.min()
        print(f"\ninbred  eff_n: {inb.min():.2f}–{inb.max():.2f}")
        print(f"outcross eff_n: {out.min():.2f}–{out.max():.2f}")
        print(f"separation: {'CLEAN' if sep else 'OVERLAP'} "
              f"(midpoint threshold eff_n = {thr:.2f})")

    # figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5),
                                   gridspec_kw=dict(width_ratios=[1, 1.4]))
    cmap = {"inbred": "#4C72B0", "outcrossed": "#DD8452"}
    d = df.reset_index(drop=True)
    ax1.bar(range(len(d)), d.eff_n_founders, color=[cmap[l] for l in d.label])
    ax1.set_xticks(range(len(d)))
    ax1.set_xticklabels([s.split("_")[0].replace("MEAJM003-", "") for s in d["sample"]],
                        rotation=45, ha="right", fontsize=8)
    ax1.axhline(1, ls=":", c="0.5", lw=1)
    ax1.set_ylabel("eff_n_founders  = 1 / Σh²")
    ax1.set_title("kMate MOI: effective # founders per individual")
    if len(inb) and len(out) and inb.max() < out.min():
        ax1.axhline(thr, ls="--", c="k", lw=1, alpha=0.6)
    handles = [plt.Rectangle((0, 0), 1, 1, color=cmap[k]) for k in cmap]
    ax1.legend(handles, [f"{k} (hapFIRE)" for k in cmap], fontsize=8)

    # full founder composition: every founder at its real frequency, coloured by
    # founder identity (same founder -> same colour across all bars). Muted
    # tab20-family palette with a SHUFFLED founder->slot mapping (fixed seed) so
    # founders exactly 60 apart don't share a colour; white borders separate the
    # meaningful (>2%) blocks.
    palette = (list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors)
               + list(plt.cm.tab20c.colors))
    perm = np.random.default_rng(0).permutation(len(founders))
    fcolor = {f: palette[perm[i] % len(palette)] for i, f in enumerate(founders)}
    for j, (_, r) in enumerate(d.iterrows()):
        h, order = comps[r["sample"]]
        bottom = 0.0
        for i in order:
            big = h[i] >= 0.02
            ax2.bar(j, h[i], bottom=bottom, width=0.85, color=fcolor[founders[i]],
                    edgecolor="white" if big else "none", linewidth=0.8 if big else 0)
            bottom += h[i]
    ax2.set_xticks(range(len(d)))
    ax2.set_xticklabels(list(d["sample"]), rotation=45, ha="right", fontsize=8)
    ax2.set_ylabel("ecotype frequency")
    ax2.set_ylim(0, 1)
    fig.tight_layout()
    png = a.out + ".png"
    fig.savefig(png, dpi=150, bbox_inches="tight")
    print(f"saved -> {png}")


if __name__ == "__main__":
    main()
