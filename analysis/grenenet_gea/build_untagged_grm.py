#!/usr/bin/env python
"""SNP-untagged non-SNP GRM(s) -- the direct test of "what does the non-SNP layer add that SNPs
can't see" (VAREXP_SELECTION_HANDOFF.md next step), instead of the confounded genome-wide
K_nonsnp (corr(K_snp,K_nonsnp)=0.998, same genealogy end to end).

Restricts K_nonsnp's own marker set (indel+SV, MAC>=5, call>=90%, from build_class_grms.py) to
markers with best_r2 < THRESH against any nearby (+/-50kb) panel SNP (build_nonsnp_tagging.py).
Two thresholds: 0.2 (the number named in the handoff) and 0.5 (secondary/more generous, since
0.2 turned out extremely strict genome-wide -- only ~0.26% of non-SNP markers survive it) so the
marker-count/power tradeoff is visible alongside the result.

Output -> analysis/grenenet_gea/varexp/untagged_grms.npz
  founders, K_untagged_r02, K_untagged_r05, n_markers (json: per-threshold total/indel/sv counts)
Env: kmate.
"""
from __future__ import annotations
import os, sys, time, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from build_class_grms import _load_chrom, _ZZt, N

OUT = f"{lib.GEA}/varexp"
CHROMS = ["chr1", "chr2", "chr3", "chr4", "chr5"]
THRESHOLDS = {"r02": 0.2, "r05": 0.5}


def main():
    t0 = time.time()
    founders = np.load(f"{lib.PROJ}/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz",
                       allow_pickle=True)["founders"].astype(str)

    Ksum = {k: np.zeros((N, N)) for k in THRESHOLDS}
    Mc = {k: 0 for k in THRESHOLDS}
    cls_counts = {k: {"indel": 0, "sv": 0} for k in THRESHOLDS}

    for cl in CHROMS:
        tag = np.load(f"{OUT}/nonsnp_tagging_{cl}.npz", allow_pickle=True)
        col_idx, cls_arr, best_r2 = tag["col_idx"], tag["cls"], tag["best_r2"]
        no_snp = ~(best_r2 >= 0)                       # NaN = no SNP in +/-50kb at all
        rl, al, vp, vc, n_alt, n_cal = _load_chrom(cl)
        for k, thresh in THRESHOLDS.items():
            mask = no_snp | (best_r2 < thresh)
            cols = col_idx[mask]
            K, m = _ZZt(vp, vc, n_alt, n_cal, cols)
            Ksum[k] += K; Mc[k] += m
            cls_counts[k]["indel"] += int((cls_arr[mask] == 1).sum())
            cls_counts[k]["sv"] += int((cls_arr[mask] == 2).sum())
        print(f"{cl}: " + ", ".join(f"{k}={Mc[k]:,}" for k in THRESHOLDS) +
              f" ({time.time()-t0:.0f}s)", flush=True)

    out = {"founders": founders}
    nmark = {}
    for k in THRESHOLDS:
        out[f"K_untagged_{k}"] = Ksum[k] / max(Mc[k], 1)
        nmark[k] = {"total": Mc[k], **cls_counts[k]}
    np.savez(f"{OUT}/untagged_grms.npz", n_markers=json.dumps(nmark), **out)
    print("n_markers:", json.dumps(nmark))
    print(f"wrote {OUT}/untagged_grms.npz ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
