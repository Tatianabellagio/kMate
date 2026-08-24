#!/usr/bin/env python
"""Build per-GENERATION per-sample AF matrices (snp / nonsnp / smallindel), float16.

Reproduces analysis/grenenet_selection/gen_matrices/ DIRECTLY from the on-cluster
af_store (no TSV re-read, no EM, no 1.6 TB rebuild). The 2,168 cohort samples
are individual collection TIMEPOINTS spanning 3 generations; here each sample's
compact per-sample AF vector is stacked into a per-generation matrix.

Validated byte-for-byte against the Jun-2026 Drive copy:
  dtype  float16   (NaN = missing)
  rows   pool_table() samples grouped by generation: gen1=897, gen2=763, gen3=508
  kinds  snp        6,237,063 records (index_snp.npz)        all SNPs
         nonsnp     2,252,583 records (index_nonsnp.npz)     all non-SNP
         smallindel 2,026,115 records (index_smallindel.npz) non-SNP, |dlen| <= 50

Outputs (--out, default analysis/grenenet_selection/gen_matrices):
  gen{g}_{kind}_af.npy   float16 [n_samples_g x n_rec_kind]
  gen{g}.rowmeta.csv     sampleid,site,plot,generation,flowerscollected,coverage
  index_smallindel.npz   chrom,pos,ref_len,alt_len for the small-indel subset (once)

Each (gen, kind) build is independent and idempotent (skips if the output exists
with the right shape), so it shards cleanly across a SLURM array and resumes
after preemption.

Usage:
  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY analysis/grenenet_selection/common/build_gen_matrices.py --gens 1 --kinds snp
  $PY analysis/grenenet_selection/common/build_gen_matrices.py            # all gens, all kinds
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

# SV = |alt_len - ref_len| > SV_MIN_BP (lib.is_sv default). smallindel = the rest
# of the non-SNP records (length change <= SV_MIN_BP).
SV_MIN_BP = 50
KINDS = ("snp", "nonsnp", "smallindel")
ROWMETA_COLS = ["sampleid", "site", "plot", "generation",
                "flowerscollected", "coverage"]
# which per-sample af_store subdir each kind reads from (smallindel is a column
# subset of the non-SNP vectors)
_SRC = {"snp": "af_snp", "nonsnp": "af_nonsnp", "smallindel": "af_nonsnp"}


def smallindel_mask(store: str):
    """Boolean selector over the non-SNP records: True for small indels."""
    ix = np.load(f"{store}/index_nonsnp.npz")
    dlen = np.abs(ix["alt_len"].astype("int64") - ix["ref_len"].astype("int64"))
    return (dlen <= SV_MIN_BP), ix


def write_smallindel_index(store: str, out: str):
    """Write index_smallindel.npz once (atomic). Returns the non-SNP subset mask."""
    mask, ix = smallindel_mask(store)
    os.makedirs(out, exist_ok=True)
    path = f"{out}/index_smallindel.npz"
    if not os.path.exists(path):
        # PID-unique tmp so concurrent array tasks never share/steal a tmp file;
        # the atomic replace is harmless to repeat (identical content).
        tmp = f"{path}.{os.getpid()}.tmp.npz"
        np.savez(tmp, chrom=ix["chrom"][mask], pos=ix["pos"][mask],
                 ref_len=ix["ref_len"][mask], alt_len=ix["alt_len"][mask])
        os.replace(tmp, path)
        print(f"  wrote {path} ({int(mask.sum()):,} records)", flush=True)
    return mask


def n_records(kind: str, store: str, simask) -> int:
    meta = json.load(open(f"{store}/meta.json"))
    if kind == "snp":
        return meta["n_snp"]
    if kind == "nonsnp":
        return meta["n_nonsnp"]
    return int(simask.sum())


def _done(path: str, n_rows: int, n_rec: int) -> bool:
    if not os.path.exists(path):
        return False
    try:
        shp = np.load(path, mmap_mode="r").shape
        return shp == (n_rows, n_rec)
    except Exception:
        return False


def build(gen: int, kind: str, store: str, out: str, simask):
    pt = lib.pool_table(store)
    sub = pt[pt.generation == gen].reset_index(drop=True)
    os.makedirs(out, exist_ok=True)
    # rowmeta (identical for all kinds of a generation) — atomic + PID-unique tmp
    # so the 3 same-gen array tasks don't corrupt it racing each other.
    rm = f"{out}/gen{gen}.rowmeta.csv"
    rm_tmp = f"{rm}.{os.getpid()}.tmp"
    sub[ROWMETA_COLS].to_csv(rm_tmp, index=False)
    os.replace(rm_tmp, rm)

    n_rec = n_records(kind, store, simask)
    out_path = f"{out}/gen{gen}_{kind}_af.npy"
    if _done(out_path, len(sub), n_rec):
        print(f"  gen{gen} {kind}: already done ({len(sub)} x {n_rec:,}) — skip",
              flush=True)
        return

    src = _SRC[kind]
    tmp = out_path + ".tmp.npy"
    mat = np.lib.format.open_memmap(tmp, mode="w+", dtype=np.float16,
                                    shape=(len(sub), n_rec))
    for i, s in enumerate(sub.sampleid):
        f = lib.decode_af(np.load(f"{store}/{src}/{s}.npy"))   # float32 + NaN
        if kind == "smallindel":
            f = f[simask]
        mat[i] = f.astype(np.float16)
        if (i + 1) % 100 == 0:
            print(f"    gen{gen} {kind}: {i+1}/{len(sub)}", flush=True)
    mat.flush()
    del mat
    os.replace(tmp, out_path)
    gb = len(sub) * n_rec * 2 / 1e9
    print(f"  gen{gen} {kind}: {len(sub)} x {n_rec:,} ({gb:.1f} GB) -> {out_path}",
          flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default=lib.AF_STORE)
    ap.add_argument("--out", default=f"{lib.GEA}/gen_matrices")
    ap.add_argument("--gens", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--kinds", nargs="+", default=list(KINDS), choices=KINDS)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    # smallindel index + mask only when a smallindel kind is requested, so snp /
    # nonsnp tasks don't all pile onto the shared index write.
    simask = (write_smallindel_index(args.store, args.out)
              if "smallindel" in args.kinds else None)
    for g in args.gens:
        for k in args.kinds:
            build(g, k, args.store, args.out, simask)


if __name__ == "__main__":
    main()
