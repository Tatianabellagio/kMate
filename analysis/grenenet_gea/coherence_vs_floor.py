#!/usr/bin/env python
"""Coherence cost of coarsening: per-window PC1-VE on the gen9 pools, for the base CLQ0.9
map and each merge-floor map (floor 8/15/25/40). Merging across LD boundaries should lower
PC1-VE; this quantifies the coherence we trade for trackability.

Output: results/grenenet_gea/blocks_mcf90/coarse/coherence_vs_floor.csv (floor,chrom,start,end,n_gen9,pc1_ve)
Run in kmate env.
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eval_block_coherence import load_allclass_af, pc1_var_explained

BR = "results/grenenet_gea/blocks_mcf90"
COARSE = f"{BR}/coarse"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
FLOORS = [0, 8, 15, 25, 40]   # 0 = base CLQ0.9 map


def win_file(chrom, floor):
    cl = chrom.lower()
    return (f"{BR}/{cl}_clq0.9_blocks_clq0.9.tsv" if floor == 0
            else f"{COARSE}/{cl}_floor{floor}.tsv")


def main():
    rows = []
    for chrom in CHROMS:
        print(f"[{chrom}] loading gen9 pools ...", flush=True)
        pos, AF = load_allclass_af(chrom, "gen9")
        for floor in FLOORS:
            w = pd.read_csv(win_file(chrom, floor), sep="\t")
            n = 0
            for s, e in zip(w.start_pos, w.end_pos):
                lo = np.searchsorted(pos, s); hi = np.searchsorted(pos, e, side="right")
                if hi - lo >= 2:
                    rows.append((floor, chrom, int(s), int(e), int(hi - lo),
                                 round(pc1_var_explained(AF[:, lo:hi]), 4)))
                    n += 1
            print(f"  floor {floor}: {n} windows scored", flush=True)
    df = pd.DataFrame(rows, columns=["floor", "chrom", "start", "end", "n_gen9", "pc1_ve"])
    out = f"{COARSE}/coherence_vs_floor.csv"
    df.to_csv(out, index=False)
    print(f"\nwrote {len(df)} windows -> {out}")
    print(f"\n{'floor':>6}{'#win':>9}{'medVar':>8}{'medVE':>8}{'%VE>=0.7':>10}{'%VE>=0.8':>10}")
    for floor in FLOORS:
        d = df[df.floor == floor]
        lab = "base" if floor == 0 else str(floor)
        print(f"{lab:>6}{len(d):>9}{int(d.n_gen9.median()):>8}{d.pc1_ve.median():>8.3f}"
              f"{100*(d.pc1_ve>=0.7).mean():>9.0f}%{100*(d.pc1_ve>=0.8).mean():>9.0f}%")


if __name__ == "__main__":
    main()
