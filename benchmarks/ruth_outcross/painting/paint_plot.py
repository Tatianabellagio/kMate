#!/usr/bin/env python3
"""Chromosome painting from kMate window-mode h_blocks.

Each window's local founder mixture h is the local ancestry of the individual.
We paint position (x) vs stacked local founder fraction (y), coloured by founder.
Inbred -> one founder flat across the chromosome. Recombinant outcross -> the
dominant founder(s) switch along the chromosome; each switch is a recombination
breakpoint. (kMate h is unphased diploid dosage: ~0.5 = heterozygous, ~1 = homozygous.)

Usage:
  paint_plot.py --tag win10k   [--chrom Chr1] [--min-show 0.12]
  paint_plot.py --tag dynld
"""
import argparse, glob, os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = Path("/global/scratch/users/tbellg/kmate")
OUT = HERE / "out"
KMER_PA_META = ROOT / "data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1.meta.npz"

FOUNDERS = np.asarray(np.load(KMER_PA_META, allow_pickle=True)["founders"]).astype(str)
# same consistent colour scheme as the composition figure (tab20 shuffled, seed 0)
PALETTE = (list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors)
           + list(plt.cm.tab20c.colors))
_perm = np.random.default_rng(0).permutation(len(FOUNDERS))
FCOLOR = {i: PALETTE[_perm[i] % len(PALETTE)] for i in range(len(FOUNDERS))}


def load_blocks(sample, tag, chrom):
    p = OUT / f"{sample}_{chrom}_{tag}.h_blocks_per_chrom.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=True)
    hb = np.asarray(z[f"{chrom}_h_blocks"], float)          # (n_blocks, n_founders)
    start = np.asarray(z[f"{chrom}_block_start"], float)
    end = np.asarray(z[f"{chrom}_block_end"], float)
    # normalise each block to a composition
    s = hb.sum(1, keepdims=True); s[s == 0] = 1
    return hb / s, start, end


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="win10k", help="win10k | dynld")
    ap.add_argument("--chrom", default="Chr1")
    ap.add_argument("--min-show", type=float, default=0.12,
                    help="founders reaching this fraction in any window get their own "
                         "colour; the rest are lumped grey")
    ap.add_argument("--manifest", default=str(HERE / "samples.tsv"))
    a = ap.parse_args()

    man = pd.read_csv(a.manifest, sep="\t").sort_values("label")
    data = {}
    for _, r in man.iterrows():
        d = load_blocks(r["sample_id"], a.tag, a.chrom)
        if d is not None:
            data[r["sample_id"]] = (d, r["label"])
    if not data:
        print(f"No h_blocks for tag={a.tag} chrom={a.chrom} yet.")
        return

    # founders that matter anywhere (across all samples) -> coloured; rest grey
    big = set()
    for (hb, _, _), _ in data.values():
        big |= set(np.where(hb.max(0) >= a.min_show)[0].tolist())
    big = sorted(big)
    print(f"[{a.tag}/{a.chrom}] coloured founders ({len(big)}): "
          f"{[FOUNDERS[i] for i in big]}")

    n = len(data)
    fig, axes = plt.subplots(n, 1, figsize=(13, 1.9 * n + 0.6), squeeze=False, sharex=True)
    for ax, (sample, ((hb, start, end), label)) in zip(axes[:, 0], data.items()):
        mb = (start + end) / 2 / 1e6
        widths = (end - start) / 1e6
        bottom = np.zeros(len(hb))
        # coloured founders first (stable order), then the grey remainder
        for i in big:
            ax.bar(mb, hb[:, i], width=widths, bottom=bottom, align="center",
                   color=FCOLOR[i], linewidth=0)
            bottom += hb[:, i]
        ax.bar(mb, 1 - bottom, width=widths, bottom=bottom, align="center",
               color="0.85", linewidth=0)
        ax.set_ylim(0, 1); ax.set_ylabel(f"{sample.replace('MEAJM003-','')}\n{label}",
                                         fontsize=8, rotation=0, ha="right", va="center")
        ax.set_yticks([0, 0.5, 1])
        ax.margins(x=0)
    axes[-1, 0].set_xlabel(f"{a.chrom} position (Mb)")
    # legend of coloured founders
    handles = [plt.Rectangle((0, 0), 1, 1, color=FCOLOR[i]) for i in big] + \
              [plt.Rectangle((0, 0), 1, 1, color="0.85")]
    fig.legend(handles, [FOUNDERS[i] for i in big] + ["other"],
               loc="upper center", ncol=min(len(big) + 1, 10), fontsize=7,
               frameon=False, bbox_to_anchor=(0.5, 1.0), title="founder (local ancestry)")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = HERE / f"painting_{a.tag}_{a.chrom}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()
