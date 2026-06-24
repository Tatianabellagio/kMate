#!/usr/bin/env python3
"""Genome-wide chromosome painting (all 5 chroms) per individual, from kMate
window-mode h_blocks, annotated with hapFIRE's call.

One figure per sample: 5 chromosome tracks (to scale), position vs stacked local
founder fraction. Expectation by hapFIRE category:
  inbred  -> one founder flat on every chromosome
  F1      -> SAME two founders at ~0.5 each across the WHOLE genome, NO switches
  recomb  -> dominant founder(s) switch along chromosomes (recombination breakpoints)
  complex -> several founders

Usage: paint_genome.py [--tag full_win10k] [--min-show 0.12]
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
OUT = HERE / "out"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
CHRLEN = {"Chr1": 30.43, "Chr2": 19.70, "Chr3": 23.46, "Chr4": 18.59, "Chr5": 26.98}

FOUNDERS = np.asarray(np.load(ROOT / "data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1.meta.npz",
                              allow_pickle=True)["founders"]).astype(str)
PALETTE = (list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors) + list(plt.cm.tab20c.colors))
_perm = np.random.default_rng(0).permutation(len(FOUNDERS))
FCOLOR = {i: PALETTE[_perm[i] % len(PALETTE)] for i in range(len(FOUNDERS))}
FID = {f: i for i, f in enumerate(FOUNDERS)}


def load(sample, tag, chrom):
    p = OUT / f"{sample}_{chrom}_{tag}.h_blocks_per_chrom.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=True)
    hb = np.asarray(z[f"{chrom}_h_blocks"], float)
    s = hb.sum(1, keepdims=True); s[s == 0] = 1
    return hb / s, np.asarray(z[f"{chrom}_block_start"], float), np.asarray(z[f"{chrom}_block_end"], float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="full_win10k")
    ap.add_argument("--min-show", type=float, default=0.12)
    a = ap.parse_args()
    man = pd.read_csv(HERE / "extremes_samples.tsv", sep="\t")

    for _, r in man.iterrows():
        sample, cat = r["sample_id"], r["category"]
        hf_ids = [x for x in str(r["hapfire_ecotypes"]).split("|") if x.isdigit()]
        hf = [FID[x] for x in hf_ids if x in FID]   # hapFIRE ecotypes as founder INDICES
        chroms = {c: load(sample, a.tag, c) for c in CHROMS}
        chroms = {c: v for c, v in chroms.items() if v is not None}
        if not chroms:
            print(f"  [skip] {sample}: no painting yet"); continue

        # founders to colour: those reaching min_show anywhere, plus hapFIRE's calls
        big = set(hf)
        for hb, _, _ in chroms.values():
            big |= set(np.where(hb.max(0) >= a.min_show)[0].tolist())
        big = sorted(big)

        fig, axes = plt.subplots(len(chroms), 1, figsize=(12, 1.25 * len(chroms) + 1.2),
                                 squeeze=False)
        xmax = max(CHRLEN.values())
        for ax, (c, (hb, st, en)) in zip(axes[:, 0], chroms.items()):
            order = np.argsort(st)                      # ensure ascending position
            mb = ((st + en) / 2 / 1e6)[order]
            # stacked filled areas (PolyCollection) — fast vs per-window bars
            layers = [hb[order, i] for i in big]
            other = 1 - np.sum([hb[order, i] for i in big], axis=0)
            layers.append(np.clip(other, 0, 1))
            cols = [FCOLOR[i] for i in big] + [(0.85, 0.85, 0.85)]
            polys = ax.stackplot(mb, layers, colors=cols, linewidth=0)
            for p in polys:
                p.set_rasterized(True)
            ax.set_xlim(0, xmax); ax.set_ylim(0, 1); ax.set_yticks([])
            ax.set_ylabel(c, fontsize=8, rotation=0, ha="right", va="center")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        axes[-1, 0].set_xlabel("position (Mb)")
        fig.suptitle(f"{sample}", fontsize=11, y=0.99)
        handles = [plt.Rectangle((0, 0), 1, 1, color=FCOLOR[i]) for i in big] + \
                  [plt.Rectangle((0, 0), 1, 1, color="0.85")]
        labels = [f"{FOUNDERS[i]}{'*' if i in set(hf) else ''}" for i in big] + ["other"]
        fig.legend(handles, labels, loc="upper center", ncol=min(len(big) + 1, 12),
                   fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.96))
        fig.tight_layout(rect=[0, 0, 1, 0.90])
        outp = HERE / f"genome_{cat}_{sample}.png"
        fig.savefig(outp, dpi=150, bbox_inches="tight")
        print(f"saved -> {outp}  (coloured founders: {[FOUNDERS[i] for i in big]})")


if __name__ == "__main__":
    main()
