#!/usr/bin/env python
"""Build the gen0 (SEEDMIX founding) counterpart of build_class_matrices.py.

Same snp/sv/smallindel split and MAF filter as the evolved-generation class
matrices, but the "pool matrix" here is the 8 raw SEEDMIX rep TSVs (founding
p0) instead of a site_gen_plot pool matrix — so gen1/gen2/gen3/gen9 (evolved)
and gen0 (founding) end up in the exact same {class}_gen{g}_af.npy /
{class}_gen{g}.records.csv layout.

Outputs (--out, default analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/class_matrices):
  {class}_gen0_af.npy       float32 [8 reps x n_kept]   (one row per SEEDMIX rep)
  {class}_gen0.records.csv  chrom,pos,ref_len,alt_len,p_bar,maf,n_finite,col
  gen0.pools.csv            rep id (SEEDMIX_S1..S8)
"""
from __future__ import annotations
import argparse, glob, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

SV_MIN_BP = 50
CLASSES = ("snp", "sv", "smallindel")


def class_mask(idx: dict, cls: str) -> np.ndarray:
    if cls == "snp":
        return np.ones(len(idx["pos"]), dtype=bool)
    dlen = np.abs(idx["alt_len"].astype("int64") - idx["ref_len"].astype("int64"))
    return dlen > SV_MIN_BP if cls == "sv" else dlen <= SV_MIN_BP


def load_seedmix_matrix(store: str, snp_mask: np.ndarray, src: str) -> np.ndarray:
    """[8 reps x n_src] AF matrix for the SNP or non-SNP subset, in index_{src}.npz order."""
    fs = sorted(f for f in glob.glob(f"{lib.SEEDMIX}/SEEDMIX_S*.tsv") if "_Chr" not in f)
    if len(fs) != 8:
        raise FileNotFoundError(f"expected 8 SEEDMIX rep TSVs, found {len(fs)} in {lib.SEEDMIX}")
    mask = snp_mask if src == "snp" else ~snp_mask
    rows = []
    for f in fs:
        af = pd.read_csv(f, sep="\t", usecols=["alt_freq"]).alt_freq.to_numpy(dtype=np.float32)
        rows.append(af[mask])
        print(f"    read {os.path.basename(f)}", flush=True)
    return np.vstack(rows), [os.path.basename(f).replace(".tsv", "") for f in fs]


def build(cls: str, store: str, out: str, maf_min: float):
    src = "snp" if cls == "snp" else "nonsnp"
    idx = dict(np.load(f"{store}/index_{src}.npz", allow_pickle=True))
    snp_mask = np.load(f"{store}/snp_mask.npy")
    mat, reps = load_seedmix_matrix(store, snp_mask, src)
    cmask = class_mask(idx, cls)
    sub = np.ascontiguousarray(mat[:, cmask])
    cls_cols = np.where(cmask)[0]

    finite = np.isfinite(sub)
    n_finite = finite.sum(axis=0)
    with np.errstate(invalid="ignore"):
        p_bar = np.nansum(np.where(finite, sub, 0.0), axis=0) / np.maximum(n_finite, 1)
    maf = np.minimum(p_bar, 1.0 - p_bar)
    keep = maf >= maf_min
    kept = np.where(keep)[0]
    print(f"  gen0 {cls}: {len(cls_cols):,} {cls} records -> {len(kept):,} kept "
          f"(MAF>={maf_min}; n_finite in [{int(n_finite.min())},{int(n_finite.max())}]/8)", flush=True)
    if len(kept) == 0:
        return

    os.makedirs(out, exist_ok=True)
    m = sub[:, kept].astype(np.float32)
    np.save(f"{out}/{cls}_gen0_af.npy", m)

    chrom = idx["chrom"][cls_cols][kept]
    pos = idx["pos"][cls_cols][kept]
    block = lib.assign_ld_blocks(np.asarray(chrom, dtype=str), np.asarray(pos))
    recs = pd.DataFrame(dict(
        chrom=chrom, pos=pos,
        ref_len=idx["ref_len"][cls_cols][kept], alt_len=idx["alt_len"][cls_cols][kept],
        p_bar=p_bar[kept], maf=maf[kept], n_finite=n_finite[kept],
        col=cls_cols[kept], block=block))
    recs.to_csv(f"{out}/{cls}_gen0.records.csv", index=False)
    pd.DataFrame({"pool": reps}).to_csv(f"{out}/gen0.pools.csv", index=False)
    print(f"    -> {out}/{cls}_gen0_af.npy [{m.shape[0]} x {m.shape[1]:,}], "
          f"{recs.block.ne('').sum():,} with LD block", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default=lib.AF_STORE)
    ap.add_argument("--out", default=f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/class_matrices")
    ap.add_argument("--classes", nargs="+", default=list(CLASSES), choices=CLASSES)
    ap.add_argument("--maf-min", type=float, default=0.05)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for c in args.classes:
        build(c, args.store, args.out, args.maf_min)


if __name__ == "__main__":
    main()
