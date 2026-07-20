#!/usr/bin/env python
"""Build filtered per-CLASS, per-GENERATION pool matrices for the phase-1 GEA.

Takes the flower-weighted `site_gen_plot` pool matrices (the phase-1 analysis
unit) and splits/filters them into the three variant classes the kMate panel
adds over the SNP-only phase-1 paper:

    snp         ref_len==1 & alt_len==1                 (from pool_gen{g}_snp)
    sv          |alt_len-ref_len| >  SV_MIN_BP (50)     (from pool_gen{g}_nonsnp)
    smallindel  non-SNP & |alt_len-ref_len| <= SV_MIN_BP(from pool_gen{g}_nonsnp)

Filters (mirroring phase-1 `maf05 mincount05`; both CLI-tunable):
  * MAF >= --maf-min (0.05)         contemporary MAF = min(p_bar, 1-p_bar),
                                    p_bar = NaN-aware mean AF across this gen's pools
  * record finite (non-NaN) in >= --min-finite-frac of the pools (panel/coverage
    availability; the kMate analog of phase-1's min-count)
  * invariant records (p_bar in {0,1} after filtering) are dropped (MAF=0 < maf-min)

Each kept record is assigned its phase-1 hapFIRE LD block (lib.assign_ld_blocks).

Outputs (--out, default analysis/grenenet_gea/phase1_replication/results/class_matrices):
  {class}_gen{g}_af.npy       float32 [n_pools_g x n_kept]   (NaN where no data)
  {class}_gen{g}.records.csv  chrom,pos,ref_len,alt_len,p_bar,maf,n_finite,
                              col (orig column index in the class),block
                              (row-aligned to the matrix columns)
  gen{g}.pools.csv            pool,site,plot,bio1,...  (row-aligned to matrix rows)

Usage:
  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY build_class_matrices.py --classes sv smallindel --gens 1 3
  $PY build_class_matrices.py --classes snp           --gens 1 3   # after snp pools build
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

SV_MIN_BP = 50
CLASSES = ("snp", "sv", "smallindel", "nonsnp")
# which pool-matrix kind + af_store index each class is carved from
_SRC = {"snp": "snp", "sv": "nonsnp", "smallindel": "nonsnp", "nonsnp": "nonsnp"}


def class_mask(kind_idx: dict, cls: str) -> np.ndarray:
    """Boolean over the records of the source pool matrix selecting class `cls`."""
    if cls in ("snp", "nonsnp"):
        # snp source is all-SNP, nonsnp source is all-non-SNP (indel+SV pooled)
        return np.ones(len(kind_idx["pos"]), dtype=bool)
    dlen = np.abs(kind_idx["alt_len"].astype("int64") - kind_idx["ref_len"].astype("int64"))
    return dlen > SV_MIN_BP if cls == "sv" else dlen <= SV_MIN_BP


def build(gen: int, cls: str, pooldir: str, store: str, out: str,
          maf_min: float, min_finite_frac: float):
    src = _SRC[cls]
    mat_path = f"{pooldir}/pool_gen{gen}_{src}_af.npy"
    meta_path = f"{pooldir}/pool_gen{gen}_{src}.meta.csv"
    if not os.path.exists(mat_path):
        print(f"  SKIP gen{gen} {cls}: missing {mat_path}", flush=True)
        return
    pmeta = pd.read_csv(meta_path)                      # pools x (pool,site,plot,bio1..)
    idx = dict(np.load(f"{store}/index_{src}.npz", allow_pickle=True))
    cmask = class_mask(idx, cls)                        # columns of source belonging to cls

    # Load the WHOLE matrix into RAM with ONE sequential read, then subset columns
    # in-memory. (A memmap + fancy column-index would do millions of strided reads
    # across the file — pathological on shared scratch. gen1 nonsnp ~2.9 GB fits.)
    af = np.load(mat_path)                              # [pools x n_src], float32, NaN=missing
    n_pools = af.shape[0]
    sub = np.ascontiguousarray(af[:, cmask])            # [pools x n_cls]  in-RAM subset
    del af
    cls_cols = np.where(cmask)[0]                       # orig column index within the source

    finite = np.isfinite(sub)
    n_finite = finite.sum(axis=0)                       # per record
    with np.errstate(invalid="ignore"):
        p_bar = np.nansum(np.where(finite, sub, 0.0), axis=0) / np.maximum(n_finite, 1)
    maf = np.minimum(p_bar, 1.0 - p_bar)

    # The finite/coverage filter is INERT on kMate data: kMate projects AF from the
    # per-chrom founder-haplotype reconstruction (h), so every record is finite in
    # every pool (n_finite == n_pools). Verified across all gens x classes (~8.5M
    # records each): min n_finite == n_pools, 0 records below the 50% threshold, so
    # MAF is the sole gatekeeper. Term commented out (kept, not deleted, so it can be
    # re-enabled via --min-finite-frac if this is ever reused on data WITH missingness,
    # e.g. short-read pools). n_finite is still computed and written to records.csv.
    keep = (maf >= maf_min)  # & (n_finite >= int(np.ceil(min_finite_frac * n_pools)))
    kept = np.where(keep)[0]
    print(f"  gen{gen} {cls}: {len(cls_cols):,} {cls} records -> "
          f"{len(kept):,} kept (MAF>={maf_min}; finite filter inert on kMate: "
          f"n_finite in [{int(n_finite.min())},{int(n_finite.max())}]/{n_pools})", flush=True)
    if len(kept) == 0:
        return

    os.makedirs(out, exist_ok=True)
    mat = sub[:, kept].astype(np.float32)
    np.save(f"{out}/{cls}_gen{gen}_af.npy", mat)

    chrom = idx["chrom"][cls_cols][kept]
    pos = idx["pos"][cls_cols][kept]
    block = lib.assign_ld_blocks(np.asarray(chrom, dtype=str), np.asarray(pos))
    recs = pd.DataFrame(dict(
        chrom=chrom, pos=pos,
        ref_len=idx["ref_len"][cls_cols][kept], alt_len=idx["alt_len"][cls_cols][kept],
        p_bar=p_bar[kept], maf=maf[kept], n_finite=n_finite[kept],
        col=cls_cols[kept], block=block))
    recs.to_csv(f"{out}/{cls}_gen{gen}.records.csv", index=False)
    # pool rowmeta (same for every class of a generation; write once)
    pmeta.to_csv(f"{out}/gen{gen}.pools.csv", index=False)
    print(f"    -> {out}/{cls}_gen{gen}_af.npy [{mat.shape[0]} x {mat.shape[1]:,}], "
          f"{recs.block.ne('').sum():,} with LD block", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pooldir", default=f"{lib.GEA}/pool_matrices")
    ap.add_argument("--store", default=lib.AF_STORE)
    ap.add_argument("--out", default=f"{lib.GEA}/phase1_replication/results/class_matrices")
    ap.add_argument("--classes", nargs="+", default=list(CLASSES), choices=CLASSES)
    ap.add_argument("--gens", type=int, nargs="+", default=[1, 3])
    ap.add_argument("--maf-min", type=float, default=0.05)
    ap.add_argument("--min-finite-frac", type=float, default=0.5)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for g in args.gens:
        for c in args.classes:
            build(g, c, args.pooldir, args.store, args.out,
                  args.maf_min, args.min_finite_frac)


if __name__ == "__main__":
    main()
