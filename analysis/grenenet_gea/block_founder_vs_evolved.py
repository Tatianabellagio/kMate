#!/usr/bin/env python
"""Per-block founder vs evolved PC1-VE (recombination diagnostic) on a given block map.
founder_ve = PC1-VE of the 231-founder genotype matrix per block; evolved_ve = PC1-VE of
the gen9 pool AF; eff_dim = participation ratio of the evolved covariance eigenvalues.
evolved_ve << founder_ve => LD broke during evolution (recombination / fine-mapping lead).

Env: BLOCKS_DIR + MINCF align to the block map. Run in kmate env.
Usage: block_founder_vs_evolved.py --chrom Chr1 --clqcut 0.9
"""
import os, sys, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from recompute_blocks import build_common_matrix
from eval_block_coherence import load_allclass_af, pc1_var_explained

BR = os.environ.get("BLOCKS_DIR", "results/grenenet_gea/blocks_recompute")
MAF = 0.05
MINCF = float(os.environ.get("MINCF", "0.5"))


def eff_dim(M):
    M = np.array(M, float)
    cm = np.nanmean(M, 0); inds = np.where(np.isnan(M)); M[inds] = np.take(cm, inds[1])
    M = M - M.mean(0); sd = M.std(0); sd[sd == 0] = 1.0; M = M / sd
    if M.shape[1] < 2:
        return np.nan
    ev = np.linalg.svd(M, compute_uv=False) ** 2
    return float((ev.sum() ** 2) / (ev ** 2).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="Chr1"); ap.add_argument("--clqcut", default="0.9")
    a = ap.parse_args(); chrlc = a.chrom.lower()
    _, raw, positions = build_common_matrix(f"panel/arch3/{chrlc}/var_pa_231_arch3_{chrlc}", MAF, MINCF)
    positions = np.asarray(positions); fgeno = raw
    ppos, AF = load_allclass_af(a.chrom, "gen9")
    bl = pd.read_csv(f"{BR}/{chrlc}_clq{a.clqcut}_blocks_clq{a.clqcut}.tsv", sep="\t")
    rows = []
    for s, e in zip(bl.start_pos, bl.end_pos):
        lo = int(np.searchsorted(positions, s)); hi = int(np.searchsorted(positions, e, side="right"))
        if hi - lo < 2:
            continue
        fve = pc1_var_explained(fgeno[:, lo:hi])
        keep = []
        for p in positions[lo:hi]:
            j = np.searchsorted(ppos, p)
            if j < len(ppos) and ppos[j] == p:
                keep.append(j)
        if len(keep) < 2:
            continue
        rows.append((a.chrom, int(s), int(e), hi - lo, round(fve, 4),
                     round(pc1_var_explained(AF[:, keep]), 4), round(eff_dim(AF[:, keep]), 3)))
    df = pd.DataFrame(rows, columns=["chrom", "start_pos", "end_pos", "n_variants",
                                     "founder_ve", "evolved_ve", "eff_dim"])
    out = f"{BR}/{chrlc}_clq{a.clqcut}_founder_vs_evolved_ve.csv"
    df.to_csv(out, index=False)
    print(f"wrote {len(df)} blocks -> {out}")


if __name__ == "__main__":
    main()
