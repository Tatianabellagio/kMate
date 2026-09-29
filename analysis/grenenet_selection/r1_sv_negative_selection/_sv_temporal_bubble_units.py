#!/usr/bin/env python
"""Per-bubble table behind `sv_temporal_bubble.ipynb`: allele count, exclusivity, p0, tested flag.

Split out of `_sv_temporal_bubble.py` (whose functions it imports) so the notebook stays load-only.
The exclusivity column is the one that decides whether a bubble's frequency is even defined as a
sum of its alleles -- see `exclusive_mask` there.

Env: kmate. Writes results/sv_adaptive/sv_temporal_bubble_units.npz.
"""
import os, sys
import numpy as np
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import lib

_sp = importlib.util.spec_from_file_location("tb", os.path.join(HERE, "_sv_temporal_bubble.py"))
tb = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(tb)

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"


def main():
    cols, M, li, k, bcls, ins, n_alt, n_cal = tb.build_bubbles()
    nbub = M.shape[0]
    p0a = np.load(f"{tb.STORE}/p0_nonsnp.npy").astype(np.float64)[cols]
    p0b = np.asarray(M @ p0a).ravel()
    cb, mb = tb.bubble_founder_stats(cols, li, nbub)
    excl = tb.exclusive_mask(cols, li, nbub)
    keep_b = ((cb >= tb.MIN_MAC) & (cb <= 231 - tb.MIN_MAC) & (mb >= tb.CALL_MIN * 231)
              & (p0b >= tb.MIN_P0) & (p0b <= 1 - tb.MIN_P0) & excl)
    keep_noexcl = ((cb >= tb.MIN_MAC) & (cb <= 231 - tb.MIN_MAC) & (mb >= tb.CALL_MIN * 231)
                   & (p0b >= tb.MIN_P0) & (p0b <= 1 - tb.MIN_P0))
    # the ALLELE filter used by the published analysis is position-based: founder_panel_keep passes a
    # record if ANY panel record at that position clears MAC >= 12. So an allele can enter the tested
    # set on a sibling allele's frequency while being private itself. n_alt is its OWN founder count.
    idx = np.load(f"{tb.STORE}/index_nonsnp.npz")
    keep_a = (lib.founder_panel_keep(idx["chrom"].astype("U5")[cols], idx["pos"].astype(np.int64)[cols],
                                     min_mac=tb.MIN_MAC) & (p0a >= tb.MIN_P0) & (p0a <= 1 - tb.MIN_P0))
    np.savez_compressed(f"{OUT}/sv_temporal_bubble_units.npz", k=k, excl=excl, p0b=p0b, cb=cb, mb=mb,
                        bcls=bcls.astype(str), keep_b=keep_b, keep_noexcl=keep_noexcl,
                        p0a=p0a.astype(np.float32), k_of_allele=k[li], ins=ins,
                        n_alt=n_alt, n_cal=n_cal, keep_a=keep_a)
    t = keep_a
    print(f"  tested ALLELES {t.sum():,}: own founder MAC < {tb.MIN_MAC} for {int((n_alt[t]<tb.MIN_MAC).sum()):,} "
          f"({100*(n_alt[t]<tb.MIN_MAC).mean():.1f}%), private (MAC==1) for {int((n_alt[t]==1).sum()):,} "
          f"({100*(n_alt[t]==1).mean():.1f}%)")
    print(f"  tested BUBBLES {keep_b.sum():,}: own founder MAC < {tb.MIN_MAC} for "
          f"{int((cb[keep_b]<tb.MIN_MAC).sum()):,} (0 by construction)")
    print(f"[wrote] {OUT}/sv_temporal_bubble_units.npz")
    print(f"  {cols.size:,} SV alleles -> {nbub:,} bubbles; exclusive {excl.sum():,} ({100*excl.mean():.1f}%)")
    print(f"  tested with the exclusivity requirement: {keep_b.sum():,}; without it: {keep_noexcl.sum():,}")


if __name__ == "__main__":
    main()
