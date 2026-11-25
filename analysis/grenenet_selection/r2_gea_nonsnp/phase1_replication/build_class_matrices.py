#!/usr/bin/env python
"""Per-class pool allele-frequency matrices, filtered: the GEA's record-arm input (step 1 of _submit_gea.sh).

Reads the flower-weighted site_gen_plot pool matrices (lib.POOL_MATRICES; generation 9 = each plot's last
generation, assemble_lastgen.py) and splits them into variant classes:

  snp         ref_len == alt_len == 1                       from pool_gen{g}_snp
  smallindel  non-SNP with |alt_len - ref_len| <= 50 bp     from pool_gen{g}_nonsnp
  sv          |alt_len - ref_len| > 50 bp                   from pool_gen{g}_nonsnp
  nonsnp      smallindel + sv                               from pool_gen{g}_nonsnp

Filters, per record:
  MAF >= --maf-min (0.05), with MAF = min(p_bar, 1 - p_bar) and p_bar = mean frequency over the generation's pools.
  Founder call rate >= --min-call-rate (0.50), i.e. called in >= 116 of the 231 founders. Where few founders are
  called, kMate's frequency is the frequency within the called founders: a different estimand that varies less
  across gardens and so inflates the GEA (multiaxis/callrate_calibration.py).

Writes in --out:
  {cls}_gen{g}_af.npy        float32 [pools x kept records]
  {cls}_gen{g}.records.csv   one row per matrix column: chrom,pos,ref_len,alt_len,p_bar,maf,n_finite,col,block
                             (col = the record's index in af_store/index_{snp,nonsnp}; block = hapFIRE LD block)
  gen{g}.pools.csv           one row per matrix row: pool,site,plot,...,bio1-19,pc1-3

Environment: kmate.  Usage: build_class_matrices.py --classes snp smallindel --gens 9 --out <dir>_arch3v2
"""
from __future__ import annotations
import argparse, glob, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

SV_MIN_BP = 50
N_FOUNDERS = 231
CLASSES = ("snp", "smallindel", "sv", "nonsnp")
SOURCE = {"snp": "snp", "smallindel": "nonsnp", "sv": "nonsnp", "nonsnp": "nonsnp"}   # pool matrix each class is cut from


def class_mask(idx: dict, cls: str) -> np.ndarray:
    """Which records of the source pool matrix belong to `cls`."""
    if cls in ("snp", "nonsnp"):
        return np.ones(len(idx["pos"]), dtype=bool)
    dlen = np.abs(idx["alt_len"].astype("int64") - idx["ref_len"].astype("int64"))
    return dlen > SV_MIN_BP if cls == "sv" else dlen <= SV_MIN_BP


def load_n_called(store: str) -> dict:
    """Founders called per record, {"snp": array, "nonsnp": array}, in af_store record order.

    A property of the panel, so any one sample's af_store/nc_{kind}/ file holds it. Positional: never joined on
    (chrom, pos, ref_len, alt_len), which is not unique for multi-allelic sites.
    """
    out = {}
    for kind in ("snp", "nonsnp"):
        f = sorted(glob.glob(f"{store}/nc_{kind}/*.npy"))
        if not f:
            raise SystemExit(f"no nc_{kind}/*.npy in {store}; run common/build_af_store.py")
        out[kind] = np.load(f[0]).astype(np.int32)
    return out


def add_climate_pcs(pools: pd.DataFrame) -> None:
    """pc1-3 of the 19 standardized bioclim variables over the pools, each z-scored, added as columns.

    Orientation: pc1 correlates positively with bio1; pc2 and pc3 give their largest-|loading| variable a positive
    loading (on gen 9 that is bio7 and bio10). The sign does not change an LFMM p-value, only how an effect is read.
    """
    m = pools[[f"bio{i}" for i in range(1, 20)]].to_numpy(float)
    z = (m - m.mean(0)) / m.std(0, ddof=0)
    u, s, vt = np.linalg.svd(z - z.mean(0), full_matrices=False)
    for k in range(3):
        score = u[:, k] * s[k]
        if k == 0:
            flip = np.corrcoef(score, m[:, 0])[0, 1] < 0
        else:
            flip = vt[k][np.argmax(np.abs(vt[k]))] < 0
        if flip:
            score = -score
        pools[f"pc{k + 1}"] = (score - score.mean()) / score.std(ddof=0)


def build(gen: int, cls: str, pooldir: str, store: str, out: str, maf_min: float, min_call_rate: float):
    src = SOURCE[cls]
    pools = pd.read_csv(f"{pooldir}/pool_gen{gen}_{src}.meta.csv")
    assert not pools.site.astype(int).isin(lib.EXCLUDED_SITES).any(), \
        f"pool_gen{gen}_{src} holds an excluded site {sorted(lib.EXCLUDED_SITES)}; rebuild the pool matrices"
    add_climate_pcs(pools)
    idx = dict(np.load(f"{store}/index_{src}.npz", allow_pickle=True))
    cols = np.where(class_mask(idx, cls))[0]                # the class's records, as indices into the source

    # one sequential read of the whole matrix, then an in-memory column subset (a memmap with fancy column
    # indexing does millions of strided reads, which is pathological on shared scratch)
    af = np.load(f"{pooldir}/pool_gen{gen}_{src}_af.npy")   # [pools x records], float32
    assert af.shape[0] == len(pools)
    sub = np.ascontiguousarray(af[:, cols])
    del af

    finite = np.isfinite(sub)
    n_finite = finite.sum(axis=0)                           # = number of pools for every kMate record
    p_bar = np.nansum(np.where(finite, sub, 0.0), axis=0) / np.maximum(n_finite, 1)
    maf = np.minimum(p_bar, 1.0 - p_bar)
    n_called = load_n_called(store)[src][cols]
    floor = int(np.ceil(min_call_rate * N_FOUNDERS))
    keep = np.where((maf >= maf_min) & (n_called >= floor))[0]
    print(f"gen{gen} {cls}: {len(cols):,} records -> {len(keep):,} kept "
          f"(MAF >= {maf_min}, called in >= {floor}/{N_FOUNDERS} founders)", flush=True)
    assert len(keep) > 0

    os.makedirs(out, exist_ok=True)
    np.save(f"{out}/{cls}_gen{gen}_af.npy", sub[:, keep].astype(np.float32))
    chrom, pos = idx["chrom"][cols][keep], idx["pos"][cols][keep]
    pd.DataFrame(dict(
        chrom=chrom, pos=pos, ref_len=idx["ref_len"][cols][keep], alt_len=idx["alt_len"][cols][keep],
        p_bar=p_bar[keep], maf=maf[keep], n_finite=n_finite[keep], col=cols[keep],
        block=lib.assign_ld_blocks(np.asarray(chrom, dtype=str), np.asarray(pos)),
    )).to_csv(f"{out}/{cls}_gen{gen}.records.csv", index=False)
    pools.to_csv(f"{out}/gen{gen}.pools.csv", index=False)  # the same for every class of a generation


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--classes", nargs="+", default=["snp", "smallindel"], choices=CLASSES)
    ap.add_argument("--gens", type=int, nargs="+", default=[9])
    ap.add_argument("--maf-min", type=float, default=0.05)
    ap.add_argument("--min-call-rate", type=float, default=0.50, help="founder call-rate floor (n_called / 231)")
    ap.add_argument("--pooldir", default=lib.POOL_MATRICES)
    ap.add_argument("--store", default=lib.AF_STORE)
    a = ap.parse_args()
    for d, what in ((a.pooldir, "--pooldir"), (a.store, "--store"), (a.out, "--out")):
        lib.check_version_dir(d, what)                      # never mix panel versions (KMATE_PANEL)
    for g in a.gens:
        for c in a.classes:
            build(g, c, a.pooldir, a.store, a.out, a.maf_min, a.min_call_rate)


if __name__ == "__main__":
    main()
