#!/usr/bin/env python
"""Inputs for an EVENT-level non-SNP GEA -- the multi-axis LFMM with SVs counted as events.

Why. The GEA tests every arch3 record. arch3 splits one SV into one record per assembly path,
so a single insertion appears as several "alleles" that differ by a few bp: its frequency is
split across them, alleles that are each rare never reach MAF 0.05 and are never tested, and
an SV whose alleles do reach it is tested several times over. The burden and purging work
collapses those records into truvari events (notebooks/sv_purging_fixed.ipynb). Counted that
way, 27,006 SV events clear MAF 0.05 -- against 21,805 events the record-level GEA touches --
and 5,369 of them are NEW: no single allele was common enough, their alleles together are.

What is identical to the GEA (r2_gea_nonsnp/phase1_replication/build_class_matrices.py +
build_lfmm_input.py): the pool frequency store (gen 9, 352 pools), the founding frequency p0,
Y = delta p = pool frequency - p0, MAF = min(p_bar, 1 - p_bar) over pools with MAF >= 0.05 the
only filter, the z-scored climate axis, and the LFMM (ridge, K = 16) run downstream.

What changes is only the unit, reused from r1_sv_negative_selection/results/sv_adaptive/
nonsnp_event_units.npz (_build_event_units.py): short indels (|dlen| <= 50) and MNPs stay one
unit per record; SVs (|dlen| > 50, the GEA's own SV threshold) are one unit per truvari event.
A unit's pool frequency and p0 are the SUMS over its member records (members are carried by
disjoint founders), capped at 1.

Two scans, as the GEA ran them: `nonsnp` (every unit, joint latent factors) and `sv` (SV events
only). The small-indel scan is unchanged at the unit level and is not repeated.

Check built in: a small-indel unit is a single record, so its delta p must equal the GEA's own
input column for that record exactly.

Writes event_gea/{cls}_gen9_units.csv and event_gea/{cls}_gen9_dp.npy (pools x units, float32).
env: kmate. Compute node (~20 GB RAM).
"""
import os, sys
import numpy as np, pandas as pd, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, GEA)
import lib                                                       # noqa: E402

UNITS = f"{GEA}/r1_sv_negative_selection/results/sv_adaptive/nonsnp_event_units.npz"
CMDIR = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/class_matrices"
MAF_MIN, SV_BP = 0.05, 50


def main():
    U = np.load(UNITS)
    ptr, ind, dlen, rep = U["indptr"], U["indices"], U["dlen"], U["rep"]
    nu = len(ptr) - 1; sizes = np.diff(ptr)
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    af = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy")   # pools x records
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    n_pools, n_rec = af.shape
    assert len(p0) == n_rec and ind.max() < n_rec
    print(f"{n_pools} pools x {n_rec:,} records -> {nu:,} units", flush=True)

    M = sp.csr_matrix((np.ones(len(ind)), (np.repeat(np.arange(nu), sizes), ind)), shape=(nu, n_rec))
    up0 = np.minimum(M @ p0.astype(np.float64), 1.0)
    # unit frequency per pool = sum of member frequencies; streamed over pools to bound memory
    pbar = np.zeros(nu)
    for i in range(n_pools):
        pbar += np.minimum(M @ af[i].astype(np.float64), 1.0)
    pbar /= n_pools
    maf = np.minimum(pbar, 1 - pbar)
    keep_all = maf >= MAF_MIN
    is_sv = np.abs(dlen) > SV_BP
    print(f"units with MAF >= {MAF_MIN}: {keep_all.sum():,} "
          f"({(keep_all & is_sv).sum():,} SV events, {(keep_all & ~is_sv).sum():,} short-indel/MNP units)", flush=True)

    ch = idx["chrom"].astype(str); pos = idx["pos"]; rl = idx["ref_len"]; al = idx["alt_len"]
    for cls, keep in (("nonsnp", keep_all), ("sv", keep_all & is_sv)):
        k = np.where(keep)[0]
        Mk = M[k]
        dp = np.empty((n_pools, len(k)), dtype=np.float32)
        for i in range(n_pools):
            dp[i] = np.minimum(Mk @ af[i].astype(np.float64), 1.0) - up0[k]
        r = rep[k]
        T = pd.DataFrame(dict(unit=k, chrom=ch[r], pos=pos[r], ref_len=rl[r], alt_len=al[r],
                              dlen=dlen[k], n_alleles=sizes[k], is_sv=is_sv[k],
                              p_bar=pbar[k], MAF=maf[k]))
        T["block"] = lib.assign_ld_blocks(T.chrom.to_numpy(str), T.pos.to_numpy())
        np.save(f"{HERE}/{cls}_gen9_dp.npy", dp)
        T.to_csv(f"{HERE}/{cls}_gen9_units.csv", index=False)
        print(f"  {cls}: {len(k):,} units -> {cls}_gen9_dp.npy / _units.csv", flush=True)

        if cls == "nonsnp":
            # CHECK: single-record short-indel units must reproduce the GEA's own delta p
            R = pd.read_csv(f"{CMDIR}/nonsnp_gen9.records.csv", usecols=["col"])
            A = np.load(f"{CMDIR}/nonsnp_gen9_af.npy", mmap_mode="r")
            colpos = pd.Series(np.arange(len(R)), index=R["col"].to_numpy())
            single = np.where((sizes[k] == 1) & ~is_sv[k])[0]
            rng = np.random.default_rng(0)
            pick = rng.choice(single, min(2000, len(single)), replace=False)
            recs = ind[ptr[k[pick]]]
            have = colpos.reindex(recs).notna().to_numpy()
            gea = np.asarray(A[:, colpos.reindex(recs[have]).astype(int).to_numpy()], np.float64) - p0[recs[have]][None, :]
            mine = dp[:, pick[have]].astype(np.float64)
            print(f"  CHECK short-indel units vs GEA input: {have.sum()} compared, "
                  f"max |diff| = {np.abs(gea - mine).max():.2e}", flush=True)


if __name__ == "__main__":
    main()
