#!/usr/bin/env python
"""AUDIT of the insertion-bubble spectra in `_sv_sfs.py` (user: "I see negative values, that is weird").

The negative values are the insertion-bubble minus SNP difference going below zero, i.e. insertion
bubbles looking LESS rare than SNPs. Two candidate explanations are tested here:

 (A) AGGREGATION. A bubble's ALT frequency is the union over its alleles, so it is the SUM of the
     allele frequencies (alleles of one bubble are mutually exclusive on a haplotype). A SNP is a
     single allele. So a k-allele bubble is compared against a 1-allele SNP and must look commoner
     the larger k is -- nothing to do with selection. Test: stratify bubbles by allele count k and
     recompute the fold rare-skew vs SNPs within each stratum. If the deficit appears only at k > 1
     and k = 1 bubbles match the SNP reference, (A) is the explanation.

 (B) CALLED-SET SEMANTICS. `_sv_sfs.py` takes `called` for a bubble as the UNION over its alleles
     (a founder counts as genotyped if genotyped at ANY allele). The stricter reading is the
     intersection. Union inflates the denominator and makes bubbles look rarer, so it works against
     the negatives, but it is quantified here anyway.

Also reports how many bubbles mix insertion and deletion records (those are split by class in
`_sv_sfs.py`, so their "bubble" is only the insertion part of the bubble).

Env: kmate. Read-only; prints a report, writes nothing.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
SV_BP, CALL_MIN, SNP_KEEP, SEED = 50, 0.9, 0.05, 0
KBINS = [(1, 1), (2, 2), (3, 4), (5, 9), (10, 10**6)]


def main():
    rng = np.random.default_rng(SEED)
    F = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    BM = np.load(f"{OUT}/arch3_bubble_map.npz")
    is_cac = F.is_cactus.to_numpy().astype(bool)
    sets = {"all": np.ones(len(F), bool), "pangenie": ~is_cac, "cactus": is_cac}

    rec = {s: dict(k=[], c=[], mu=[], mi=[], fa=[], snp=[]) for s in sets}
    mixed = both = 0
    for ch in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{ch.lower()}/var_pa_231_arch3_{ch.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64); dl = al - rl
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        bmap = pd.Series(BM[f"{ch}_bubble"], index=BM[f"{ch}_rec"])
        ins = np.where(dl > SV_BP)[0]; dele = np.where(-dl > SV_BP)[0]
        bi = set(bmap.reindex(ins).to_numpy().tolist()); bd = set(bmap.reindex(dele).to_numpy().tolist())
        mixed += len(bi & bd); both += len(bi | bd)
        snp_cols = np.where((rl == 1) & (al == 1) & (rng.random(dl.size) < SNP_KEEP))[0]

        b = bmap.reindex(ins).to_numpy().astype(np.int64)
        assert (b >= 0).all(), f"{ch}: unmatched bubble codes present"
        ub, li = np.unique(b, return_inverse=True); li = li.ravel()
        M = sp.csr_matrix((np.ones(ins.size), (li, np.arange(ins.size))), shape=(ub.size, ins.size))
        kcount = np.bincount(li, minlength=ub.size)

        for s, fm in sets.items():
            P = np.asarray(vp[fm][:, ins].todense()).T.astype(bool)
            C = np.asarray(vc[fm][:, ins].todense()).T.astype(bool)
            cb = np.asarray((M @ P) > 0)                       # carriers: union (as in _sv_sfs.py)
            mu = np.asarray((M @ C) > 0)                       # called: union   (as in _sv_sfs.py)
            mi = np.asarray((M @ (~C)) == 0)                   # called: intersection (stricter)
            rec[s]["k"].append(kcount); rec[s]["c"].append(cb.sum(1))
            rec[s]["mu"].append(mu.sum(1)); rec[s]["mi"].append(mi.sum(1))
            fa = np.asarray(vp[fm][:, ins].todense()).T.astype(bool).sum(1)
            ca = np.asarray(vc[fm][:, ins].todense()).T.astype(bool).sum(1)
            rec[s]["fa"].append(np.c_[fa, ca, b])              # per-ALLELE, for the k=1 cross-check
            Ps = np.asarray(vp[fm][:, snp_cols].todense()).T.astype(bool)
            Cs = np.asarray(vc[fm][:, snp_cols].todense()).T.astype(bool)
            rec[s]["snp"].append(np.c_[Ps.sum(1), Cs.sum(1)])
        print(f"[{ch}] done", flush=True)

    print(f"\nbubbles containing BOTH insertion and deletion records: {mixed:,} of {both:,} "
          f"({100*mixed/both:.1f}%) -- _sv_sfs.py splits these by class")

    for s, m in sets.items():
        n = int(m.sum())
        k = np.concatenate(rec[s]["k"]); c = np.concatenate(rec[s]["c"])
        mu = np.concatenate(rec[s]["mu"]); mi = np.concatenate(rec[s]["mi"])
        S = np.concatenate(rec[s]["snp"]); A = np.concatenate(rec[s]["fa"])
        keep_s = (S[:, 1] >= CALL_MIN * n) & (S[:, 0] > 0)
        fs = S[keep_s, 0] / S[keep_s, 1]; ref = (fs <= 0.05).mean()
        print(f"\n=== {s} ({n} founders) -- SNP reference: {100*ref:.1f}% at <= 5% (n={keep_s.sum():,})")
        for lab, mm in (("called = UNION (what _sv_sfs.py does)", mu), ("called = INTERSECTION (stricter)", mi)):
            keep = (mm >= CALL_MIN * n) & (c > 0)
            f = c[keep] / mm[keep]
            print(f"  {lab:<42} bubbles n={keep.sum():>7,}  <=5% {100*(f<=0.05).mean():>5.1f}%  "
                  f"fold {(f<=0.05).mean()/ref:>5.2f}x")
        print(f"  {'-- by number of alleles in the bubble (called = union) --':<42}")
        print(f"  {'k alleles':<14}{'bubbles':>9}{'<=5%':>8}{'fold vs SNP':>13}{'mean freq':>11}")
        for lo, hi in KBINS:
            sel = (k >= lo) & (k <= hi) & (mu >= CALL_MIN * n) & (c > 0)
            if sel.sum() < 50:
                continue
            f = c[sel] / mu[sel]
            lab = f"{lo}" if lo == hi else (f"{lo}-{hi}" if hi < 10**5 else f"{lo}+")
            print(f"  {lab:<14}{sel.sum():>9,}{100*(f<=0.05).mean():>7.1f}%{(f<=0.05).mean()/ref:>12.2f}x{f.mean():>11.3f}")
        ka = (A[:, 1] >= CALL_MIN * n) & (A[:, 0] > 0)
        fal = A[ka, 0] / A[ka, 1]
        print(f"  {'ALL alleles (the allele unit)':<14}{ka.sum():>9,}{100*(fal<=0.05).mean():>7.1f}%"
              f"{(fal<=0.05).mean()/ref:>12.2f}x{fal.mean():>11.3f}")


if __name__ == "__main__":
    main()
