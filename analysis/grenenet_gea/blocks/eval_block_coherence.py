#!/usr/bin/env python
"""Evaluate coherence of a recomputed block map: assign gen9 pool-AF records
(all classes) to each block by position, compute PC1 variance-explained per block
(the rank-1/single-haplotype diagnostic), and summarize vs block size.

Higher median PC1-VE = more coherent (single-haplotype) blocks = better test units.
Compare across the r2-sweep block maps to pick the coherence threshold.

Run in the kmate env. Usage:
  eval_block_coherence.py --blocks <blocks.tsv> [--blocks ...] --chrom Chr1
"""
import argparse, numpy as np, pandas as pd

CM = "analysis/grenenet_gea/r2_gea_nonsnp/phase1_replication/results/class_matrices"
CLASSES = ["snp", "smallindel", "sv"]


def load_allclass_af(chrom, gen="gen9"):
    """Return (pos sorted asc, AF matrix [pools x records] aligned to pos)."""
    poss, afs = [], []
    for c in CLASSES:
        rec = pd.read_csv(f"{CM}/{c}_{gen}.records.csv", usecols=["chrom", "pos"])
        af = np.load(f"{CM}/{c}_{gen}_af.npy", mmap_mode="r")   # pools x records
        m = (rec["chrom"] == chrom).values
        poss.append(rec["pos"].values[m])
        afs.append(np.asarray(af[:, m]))
    pos = np.concatenate(poss)
    AF = np.concatenate(afs, axis=1)            # pools x Nrec
    order = np.argsort(pos, kind="stable")
    return pos[order], AF[:, order]


def pc1_var_explained(M):
    """M: pools x variants. NaN-impute per-variant mean, standardize, PC1 VE."""
    M = np.array(M, dtype=np.float64)
    col_mean = np.nanmean(M, axis=0)
    inds = np.where(np.isnan(M))
    M[inds] = np.take(col_mean, inds[1])
    M = M - M.mean(0)
    sd = M.std(0); sd[sd == 0] = 1.0
    M = M / sd
    if M.shape[1] < 2:
        return np.nan
    s = np.linalg.svd(M, compute_uv=False)
    return float(s[0] ** 2 / np.sum(s ** 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", action="append", required=True,
                    help="block TSV (chrom start_pos end_pos n_variants); repeatable")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--gen", default="gen9")
    a = ap.parse_args()

    print(f"loading all-class {a.gen} pool AF for {a.chrom} ...")
    pos, AF = load_allclass_af(a.chrom, a.gen)
    print(f"  {AF.shape[0]} pools x {AF.shape[1]} records")

    print(f"\n{'block map':40s} {'nblk':>6} {'medVE':>7} {'VE>=0.7':>8} {'medSNP':>7}")
    for bf in a.blocks:
        bl = pd.read_csv(bf, sep="\t")
        ves, sizes = [], []
        for s, e in zip(bl["start_pos"], bl["end_pos"]):
            lo, hi = np.searchsorted(pos, s), np.searchsorted(pos, e, side="right")
            if hi - lo >= 2:
                ves.append(pc1_var_explained(AF[:, lo:hi]))
                sizes.append(hi - lo)
        ves = np.array(ves); sizes = np.array(sizes)
        name = bf.split("/")[-1]
        print(f"{name:40s} {len(bl):>6} {np.nanmedian(ves):>7.3f} "
              f"{np.mean(ves >= 0.7):>8.2f} {int(np.median(sizes)):>7d}")


if __name__ == "__main__":
    main()
