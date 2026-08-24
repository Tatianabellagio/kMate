#!/usr/bin/env python
"""Compute the per-generation allele-frequency spectrum (SFS) at site 4, split by
variant class (SNP / small indel / SV), from seedmix (gen 0) through gen 1-3.

Site 4 = the "hot" pilot site (SITE_CLIMATE). Uses the compact af_store (per-sample
uint16-encoded AF vectors) rather than the raw per-sample TSVs, so this runs in
seconds. Reads only the ~57 site-4 evolved samples + the founding p0 already cached
in the store (mean over the 8 SEEDMIX reps).

Output: analysis/grenenet_gea/sfs_time_site4.npz
  <class>_<gen>   float32 array of per-record alt_freq, class in {snp,indel,sv},
                  gen in {p0,g1,g2,g3}. NaN dropped (record not observed at that
                  timepoint in any site-4 sample of that generation).
  <class>_<gen>_n   int, sample count contributing to that generation's mean.
Also writes sfs_time_site4_summary.csv: one row per (class, gen) with n_records,
n_samples, mean/median/quantiles.
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
import lib

GEA = lib.GEA
STORE = f"{GEA}/af_store"
SITE = 4
OUT_NPZ = f"{GEA}/sfs_time_site4.npz"
OUT_CSV = f"{GEA}/sfs_time_site4_summary.csv"


def gen_mean(sample_ids: list[str], subdir: str, n: int) -> tuple[np.ndarray, int]:
    """NaN-aware mean alt_freq across samples, for one af_store subdir (af_snp/af_nonsnp)."""
    tot = np.zeros(n, dtype=np.float64)
    cnt = np.zeros(n, dtype=np.float64)
    used = 0
    for s in sample_ids:
        p = f"{STORE}/{subdir}/{s}.npy"
        if not os.path.exists(p):
            continue
        a = lib.decode_af(np.load(p))
        ok = np.isfinite(a)
        tot[ok] += a[ok]
        cnt[ok] += 1
        used += 1
    with np.errstate(invalid="ignore"):
        mean = tot / np.where(cnt > 0, cnt, np.nan)
    return mean.astype(np.float32), used


def main():
    samples = lib.list_samples()
    m = lib.sample_map(samples)
    site = m[m.site == SITE]
    print(f"site {SITE}: {len(site)} evolved samples, generations "
          f"{sorted(site.generation.unique())}")

    idx_nonsnp = np.load(f"{STORE}/index_nonsnp.npz", allow_pickle=True)
    ridx = pd.DataFrame({k: idx_nonsnp[k] for k in idx_nonsnp.files})
    is_sv = lib.sv_size(ridx).to_numpy() > lib.SV_MIN_BP   # True=SV, False=small indel
    n_nonsnp = len(ridx)
    n_snp = len(np.load(f"{STORE}/index_snp.npz", allow_pickle=True)["pos"])
    print(f"nonsnp records: {n_nonsnp:,} (SV {is_sv.sum():,} / indel {(~is_sv).sum():,}); "
          f"snp records: {n_snp:,}")

    # gen 0 = seedmix p0, already computed correctly (mean over 8 SEEDMIX reps)
    p0_nonsnp = np.load(f"{STORE}/p0_nonsnp.npy")
    p0_snp = np.load(f"{STORE}/p0_snp.npy")

    out = {}
    rows = []

    def stash(cls: str, gen: str, vals: np.ndarray, n_samp: int):
        v = vals[np.isfinite(vals)]
        out[f"{cls}_{gen}"] = v
        out[f"{cls}_{gen}_n"] = np.array(n_samp)
        rows.append(dict(cls=cls, gen=gen, n_records=len(v), n_samples=n_samp,
                          mean=float(np.mean(v)) if len(v) else np.nan,
                          median=float(np.median(v)) if len(v) else np.nan,
                          q10=float(np.quantile(v, 0.10)) if len(v) else np.nan,
                          q90=float(np.quantile(v, 0.90)) if len(v) else np.nan,
                          frac_lost=float((v < 1e-6).mean()) if len(v) else np.nan,
                          frac_fixed=float((v > 1 - 1e-6).mean()) if len(v) else np.nan))

    stash("snp", "p0", p0_snp, 8)
    stash("indel", "p0", p0_nonsnp[~is_sv], 8)
    stash("sv", "p0", p0_nonsnp[is_sv], 8)

    for g in sorted(site.generation.unique()):
        gs = site[site.generation == g].index.tolist()
        snp_mean, n_snp_used = gen_mean(gs, "af_snp", n_snp)
        nonsnp_mean, n_nonsnp_used = gen_mean(gs, "af_nonsnp", n_nonsnp)
        assert n_snp_used == n_nonsnp_used == len(gs), (
            f"gen {g}: expected {len(gs)} samples, got snp={n_snp_used} nonsnp={n_nonsnp_used}")
        tag = f"g{g}"
        stash("snp", tag, snp_mean, len(gs))
        stash("indel", tag, nonsnp_mean[~is_sv], len(gs))
        stash("sv", tag, nonsnp_mean[is_sv], len(gs))
        print(f"gen {g}: {len(gs)} samples")

    np.savez(OUT_NPZ, **out)
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_CSV, index=False)
    print(summary.to_string(index=False))
    print(f"[saved] {OUT_NPZ}")
    print(f"[saved] {OUT_CSV}")


if __name__ == "__main__":
    main()
