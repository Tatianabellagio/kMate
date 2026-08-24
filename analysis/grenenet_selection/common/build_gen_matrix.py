#!/usr/bin/env python
"""Build per-GENERATION, per-TIMEPOINT stacked AF matrices (uint16), by variant class.

Mirrors the existing gen_matrices/gen1_nonsnp_af.npy exactly: ROWS are the
collection TIMEPOINTS of one generation (usesample rows present in the af_store),
COLUMNS are the panel records of one variant class, VALUES are the af_store's
uint16 AF copied VERBATIM (scale 10000, NaN sentinel 65535) — no decode/re-encode,
so the matrix is bit-identical to stacking the per-sample vectors.

CLASSES (--kind):
  snp         all n_snp SNP records              <- af_snp/<s>.npy  (verbatim)
  smallindel  non-SNP with |alt_len-ref_len|<=50 <- af_nonsnp/<s>.npy[:, smallindel_mask]
              (the exact complement of the SV mask |dlen|>50; INCLUDES equal-length
               MNPs, dlen==0, per the analysis decision 2026-06-09)
  nonsnp      all n_nonsnp non-SNP records        <- af_nonsnp/<s>.npy (verbatim;
              rebuilds/verifies the existing file)

ROW ORDER: sorted(sampleid) within the generation. Verified to reproduce the
existing gen1_nonsnp_af.npy. The SAME order is used for every kind so the three
class matrices of a generation are row-aligned (and align to the existing nonsnp).

OUTPUTS (--out, default results_grenenet_gea/gen_matrices):
  gen{g}_{kind}_af.npy     uint16 [n_timepoints_g x n_cols_kind]
  gen{g}.rowmeta.csv       sampleid,site,plot,generation,flowerscollected,coverage
                           (row-aligned; shared across kinds, written once per gen)
  index_smallindel.npz     chrom,pos,ref_len,alt_len for smallindel columns (once)

Usage (run with the env that wrote the store, numpy>=2, OR see note on samples.npy):
  PY=/Users/tatiana/mambaforge/envs/simulations/bin/python
  $PY code/build_gen_matrix.py --kind snp smallindel --gens 1 2 3
"""
from __future__ import annotations
import argparse, glob, json, os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                       # kmate_gea_export/
STORE = f"{ROOT}/results_grenenet_gea/af_store"
OUT = f"{ROOT}/results_grenenet_gea/gen_matrices"
SAMPLES_DATA = f"{ROOT}/external/samples_data_fix57.csv"
SV_MIN_BP = 50                                     # |dlen|>SV_MIN_BP == SV (lib.is_sv)


def store_sample_ids() -> set[str]:
    """Sample ids present in the af_store (from the per-sample vector filenames;
    samples.npy is a numpy>=2 object pickle and may be unreadable under numpy<2).
    af_snp and af_nonsnp hold the SAME sample set, so use whichever is present."""
    for sd in ("af_nonsnp", "af_snp"):
        fs = glob.glob(f"{STORE}/{sd}/*.npy")
        if fs:
            return {os.path.basename(p)[:-4] for p in fs}
    return set()


def gen_rows(gen: int) -> pd.DataFrame:
    """Row table for one generation: usesample timepoints in the store, sorted by
    sampleid. Columns: sampleid, site, plot, generation, flowerscollected, coverage."""
    sd = pd.read_csv(SAMPLES_DATA)
    sd = sd[sd["usesample"]].copy()
    sd["generation"] = sd["fix_57_generation"].astype(int)
    have = store_sample_ids()
    sd = sd[sd.sampleid.isin(have) & (sd.generation == gen)]
    cols = ["sampleid", "site", "plot", "generation", "flowerscollected", "coverage"]
    return sd[cols].sort_values("sampleid").reset_index(drop=True)


def smallindel_mask() -> np.ndarray:
    """Boolean over the n_nonsnp columns: |alt_len-ref_len| <= SV_MIN_BP (incl MNPs)."""
    idx = np.load(f"{STORE}/index_nonsnp.npz")
    dlen = np.abs(idx["alt_len"].astype(np.int64) - idx["ref_len"].astype(np.int64))
    return dlen <= SV_MIN_BP


def save_smallindel_index(out: str):
    idx = np.load(f"{STORE}/index_nonsnp.npz")
    m = smallindel_mask()
    np.savez(f"{out}/index_smallindel.npz", chrom=idx["chrom"][m], pos=idx["pos"][m],
             ref_len=idx["ref_len"][m], alt_len=idx["alt_len"][m])
    print(f"  index_smallindel.npz: {int(m.sum()):,} columns")


def build(gen: int, kind: str, out: str):
    rows = gen_rows(gen)
    ids = rows.sampleid.tolist()
    rows.to_csv(f"{out}/gen{gen}.rowmeta.csv", index=False)

    if kind == "snp":
        src_dir, col_mask = "af_snp", None
    elif kind == "nonsnp":
        src_dir, col_mask = "af_nonsnp", None
    elif kind == "smallindel":
        src_dir, col_mask = "af_nonsnp", smallindel_mask()
    else:
        raise ValueError(kind)

    n_cols = int(col_mask.sum()) if col_mask is not None else \
        np.load(f"{STORE}/{src_dir}/{ids[0]}.npy", mmap_mode="r").shape[0]
    path = f"{out}/gen{gen}_{kind}_af.npy"
    tmp = path + ".tmp.npy"
    mat = np.lib.format.open_memmap(tmp, mode="w+", dtype=np.uint16,
                                    shape=(len(ids), n_cols))
    for i, s in enumerate(ids):
        v = np.load(f"{STORE}/{src_dir}/{s}.npy")
        mat[i] = v[col_mask] if col_mask is not None else v
        if (i + 1) % 100 == 0:
            print(f"  gen{gen} {kind}: {i+1}/{len(ids)} rows", flush=True)
    mat.flush(); del mat
    os.replace(tmp, path)                          # atomic: no partial file on crash
    print(f"gen{gen} {kind}: {len(ids)} x {n_cols:,} uint16 "
          f"({len(ids)*n_cols*2/1e9:.2f} GB) -> {path}")


def main():
    global STORE, SAMPLES_DATA
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default=STORE,
                    help="af_store dir (default: Drive; pass a local stage for speed)")
    ap.add_argument("--samples-data", default=SAMPLES_DATA)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--kind", nargs="+", default=["snp", "smallindel"],
                    choices=["snp", "smallindel", "nonsnp"])
    ap.add_argument("--gens", type=int, nargs="+", default=[1, 2, 3])
    args = ap.parse_args()
    STORE, SAMPLES_DATA = args.store, args.samples_data
    os.makedirs(args.out, exist_ok=True)
    if "smallindel" in args.kind:
        save_smallindel_index(args.out)
    for g in args.gens:
        for k in args.kind:
            build(g, k, args.out)


if __name__ == "__main__":
    main()
