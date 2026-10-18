#!/usr/bin/env python
"""78-founder panels for the long-read-only benchmark (kMate vs hapFIRE, and kMate's 231 panel on the same reads).

The 80 long-read assembly founders of p80 minus the two that are absent from the 231-founder panel and from the GrENE-Net
VCF (6910, 7275), so every arm has the same founder set and every simulated pool can be built and scored for all arms.
From the p80 matrices (same records, same k-mers; nothing re-derived from reads):
  data/var_pa_p78.{var_pa,var_called,meta}.npz   founder rows subset; records no longer polymorphic among the 78
                                                 (no alt carrier, or every called founder carries the alt) dropped;
                                                 data/var_pa_p78.p80_index.npy = each kept record's row in var_pa_p80
                                                 (to subset the p80 simulation truth)
  data/kmer_pa_p78_raw/kmer_pa_Chr1.*            no private k-mer filter (k-mers carried by none of the 78 dropped:
                                                 they can never match a pooled read)
  data/kmer_pa_p78_filt2inv/kmer_pa_Chr1.*       the production filter recomputed for F = 78: keep 2 <= carriers <= 77
                                                 (src/kmate/filter_kmer_pa_production.py)
Env: kmate. Compute node.
"""
import sys
import numpy as np, scipy.sparse as sp

D = "/global/scratch/projects/fc_moilab/tbellg/kmate/benchmarks/p80/data"
sys.path.insert(0, "/global/scratch/projects/fc_moilab/tbellg/kmate/src/kmate")
from filter_kmer_pa_production import production_keep_mask
DROP = {"6910", "7275"}


def main():
    vm = np.load(f"{D}/var_pa_p80.meta.npz", allow_pickle=True)
    fo = np.array([str(x) for x in vm["founders"]]); keep_f = np.array([f not in DROP for f in fo])
    assert keep_f.sum() == 78, keep_f.sum()
    vp = sp.load_npz(f"{D}/var_pa_p80.var_pa.npz").tocsr()[keep_f]; vc = sp.load_npz(f"{D}/var_pa_p80.var_called.npz").tocsr()[keep_f]
    n_call = np.asarray((vc > 0).sum(0)).ravel(); n_alt = np.asarray(((vp > 0).multiply(vc > 0)).sum(0)).ravel()
    poly = (n_alt >= 1) & (n_alt < n_call)
    idx = np.flatnonzero(poly)
    print(f"var_pa: {vp.shape[1]:,} records -> {poly.sum():,} polymorphic among the 78 ({(~poly).sum():,} dropped)", flush=True)
    sp.save_npz(f"{D}/var_pa_p78.var_pa.npz", vp.tocsc()[:, idx].tocsr())
    sp.save_npz(f"{D}/var_pa_p78.var_called.npz", vc.tocsc()[:, idx].tocsr())
    np.savez(f"{D}/var_pa_p78.meta.npz", founders=vm["founders"][keep_f],
             **{k: vm[k][idx] for k in vm.keys() if k != "founders"})
    np.save(f"{D}/var_pa_p78.p80_index.npy", idx)
    import os
    km = np.load(f"{D}/kmer_pa_p80_ours/kmer_pa_Chr1.meta.npz", allow_pickle=True)
    assert [str(x) for x in km["founders"]] == list(fo), "k-mer and variant founder order differ"
    K = sp.load_npz(f"{D}/kmer_pa_p80_ours/kmer_pa_Chr1.kmer_pa.npz").tocsr()[keep_f].tocsc()
    ac = np.asarray((K > 0).sum(0)).ravel()
    for tag, mask in (("raw", ac >= 1), ("filt2inv", production_keep_mask(ac, 78))):
        out = f"{D}/kmer_pa_p78_{tag}"; os.makedirs(out, exist_ok=True)
        cols = np.flatnonzero(mask)
        sp.save_npz(f"{out}/kmer_pa_Chr1.kmer_pa.npz", K[:, cols].tocsr())
        np.savez(f"{out}/kmer_pa_Chr1.meta.npz", kmer_index=km["kmer_index"][cols], bubble_id=km["bubble_id"][cols],
                 bubble_chrom=km["bubble_chrom"], bubble_start=km["bubble_start"], bubble_end=km["bubble_end"],
                 founders=km["founders"][keep_f])
        print(f"kmer_pa_p78_{tag}: {K.shape[1]:,} k-mers -> {len(cols):,}", flush=True)


if __name__ == "__main__":
    main()
