#!/usr/bin/env python
"""Collapse the per-sample AF TSVs into a compact per-sample NPY store.

WHY: every one of the 2,168 cohort TSVs is the SAME fixed 8.49M-record panel in
the SAME order (369 MB text each, ~800 GB total, ~1.6 TB with per-chrom files).
The record metadata (chrom/pos/ref_len/alt_len) is byte-identical across all
samples, so it is stored ONCE; alt_freq becomes a per-sample float32 vector.

DESIGN (per the analysis plan):
  - the redundant record INDEX is written once, SPLIT into SNP vs non-SNP;
  - each sample's AF (and n_called) is written as its OWN vector, also SPLIT into
    SNP vs non-SNP files — so non-SNP (SV-GEA) work never loads the heavy SNP
    bulk (73.5% / 6.24M records). SNPs are kept (separately) but out of the way.
  - per-sample vectors are the canonical compact store. Per-GENERATION matrices
    (gen1/2/3) are built from them LATER, since the 2,168 samples span 3 gens and
    don't belong stacked in one matrix.

STORE LAYOUT (--out-dir):
  meta.json            n_full, n_snp, n_nonsnp, n_samples, base, dtypes
  index_snp.npz        chrom,pos,ref_len,alt_len for SNP rows        [once]
  index_nonsnp.npz     chrom,pos,ref_len,alt_len for non-SNP rows    [once]
  snp_mask.npy         bool over the full 8.49M panel (row selector) [once]
  samples.npy          all sample IDs present in base
  af_nonsnp/<s>.npy    uint16 [n_nonsnp]   <- SV-GEA substrate
  af_snp/<s>.npy       uint16 [n_snp]      <- kept, separate, heavy
  nc_nonsnp/<s>.npy    uint8  [n_nonsnp]   <- coverage / reliability
  nc_snp/<s>.npy       uint8  [n_snp]

AF is QUANTIZED to 4 decimals: stored as uint16 = round(alt_freq * 10000), so
af in [0,1] -> [0,10000]; the value 65535 is a NaN/missing sentinel. Decode with
  af = np.where(u == 65535, np.nan, u / 10000.0)
4 decimals is ample for AF (the TSV itself carried 5) and halves AF storage vs
float32. n_called maxes at the founder count (231) so it stores as uint8.

Each per-sample non-SNP AF vector is ~4.5 MB (2.25M u16); SNP ~12.5 MB. Whole
store ~55 GB vs ~800 GB of text finals; the per-chrom TSVs (~800 GB more) can be
dropped once this is verified.

SUBCOMMANDS:
  init    read sample[0] -> SNP/non-SNP masks, split index, sample list
  convert --start S --count C : convert samples [S, S+C) into the 4 vector dirs
                                (skips any sample whose 4 vectors already exist)
"""
from __future__ import annotations
import argparse, glob, json, os
import numpy as np
import pandas as pd

PROJ = "/global/scratch/users/tbellg/kmate"
# Repointed 2026-07-07 to the full-panel-Kf_w + haploblock (--unit chrom) rerun.
OUT = f"{PROJ}/analysis/grenenet_gea/common/rerun_kfw_hb/evolved"
N_FULL_EXPECTED = 8_489_646   # segregating-only arch3 panel (sanity guard)
SUBDIRS = ("af_nonsnp", "af_snp", "nc_nonsnp", "nc_snp")
AF_SCALE = 10000        # alt_freq stored as uint16 round(af * AF_SCALE)
AF_NAN = 65535          # uint16 sentinel for NaN/missing AF (valid range 0..10000)


def encode_af(af: np.ndarray) -> np.ndarray:
    """alt_freq float -> uint16 quantized to 4 decimals; NaN -> AF_NAN sentinel."""
    finite = np.isfinite(af)
    q = np.empty(af.shape, dtype=np.uint16)
    q[~finite] = AF_NAN
    vals = np.rint(np.clip(af[finite], 0.0, 1.0) * AF_SCALE).astype(np.uint16)
    q[finite] = vals
    return q


def decode_af(u: np.ndarray) -> np.ndarray:
    """uint16 store -> float alt_freq in [0,1]; AF_NAN sentinel -> NaN."""
    f = u.astype(np.float32) / AF_SCALE
    f[u == AF_NAN] = np.nan
    return f


def list_samples(base: str) -> list[str]:
    """Genome-wide per-sample TSVs in `base` (excludes per-chrom files), sorted."""
    return sorted(os.path.basename(p)[:-4] for p in glob.glob(f"{base}/*.tsv")
                  if "_Chr" not in os.path.basename(p))


def cmd_init(args):
    base, d = args.base, args.out_dir
    for sd in SUBDIRS:
        os.makedirs(f"{d}/{sd}", exist_ok=True)
    samples = list_samples(base)
    if not samples:
        raise SystemExit(f"no per-sample TSVs in {base}")
    m = pd.read_csv(f"{base}/{samples[0]}.tsv", sep="\t",
                    usecols=["chrom", "pos", "ref_len", "alt_len"])
    n_full = len(m)
    if n_full != N_FULL_EXPECTED:
        print(f"WARNING: sample[0] has {n_full:,} rows, expected {N_FULL_EXPECTED:,}")
    snp = ((m.ref_len == 1) & (m.alt_len == 1)).to_numpy()
    chrom = m.chrom.to_numpy().astype("U5")
    pos = m.pos.to_numpy(); rl = m.ref_len.to_numpy(); al = m.alt_len.to_numpy()
    np.savez(f"{d}/index_snp.npz", chrom=chrom[snp], pos=pos[snp],
             ref_len=rl[snp], alt_len=al[snp])
    np.savez(f"{d}/index_nonsnp.npz", chrom=chrom[~snp], pos=pos[~snp],
             ref_len=rl[~snp], alt_len=al[~snp])
    np.save(f"{d}/snp_mask.npy", snp)
    np.save(f"{d}/samples.npy", np.array(samples, dtype=object), allow_pickle=True)
    meta = dict(n_full=int(n_full), n_snp=int(snp.sum()),
                n_nonsnp=int((~snp).sum()), n_samples=len(samples),
                base=base, af_dtype="uint16", nc_dtype="uint8",
                af_scale=AF_SCALE, af_nan=AF_NAN)
    with open(f"{d}/meta.json", "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"init: {meta['n_samples']} samples | SNP {meta['n_snp']:,} | "
          f"non-SNP {meta['n_nonsnp']:,} | full {n_full:,} -> {d}")


def _load_meta(d: str):
    with open(f"{d}/meta.json") as fh:
        meta = json.load(fh)
    snp = np.load(f"{d}/snp_mask.npy")
    samples = np.load(f"{d}/samples.npy", allow_pickle=True)
    return meta, snp, samples


def _paths(d: str, name: str):
    return {sd: f"{d}/{sd}/{name}.npy" for sd in SUBDIRS}


def _done(d: str, name: str, n_snp: int, n_nonsnp: int) -> bool:
    p = _paths(d, name)
    if not all(os.path.exists(v) for v in p.values()):
        return False
    try:
        return (np.load(p["af_nonsnp"], mmap_mode="r").shape[0] == n_nonsnp and
                np.load(p["af_snp"], mmap_mode="r").shape[0] == n_snp)
    except Exception:
        return False


def cmd_convert(args):
    d = args.out_dir
    meta, snp, samples = _load_meta(d)
    base, n_full = meta["base"], meta["n_full"]
    n_snp, n_nonsnp = meta["n_snp"], meta["n_nonsnp"]
    nonsnp = ~snp
    s, c = args.start, args.count
    end = min(s + c, len(samples))
    if s >= len(samples):
        print(f"convert start {s} >= n_samples {len(samples)} — nothing to do")
        return
    for j in range(s, end):
        name = str(samples[j])
        if _done(d, name, n_snp, n_nonsnp):
            print(f"  [{j}] {name} already done — skipping", flush=True)
            continue
        t = pd.read_csv(f"{base}/{name}.tsv", sep="\t",
                        usecols=["alt_freq", "n_called"])
        if len(t) != n_full:
            raise SystemExit(f"{name}: {len(t):,} rows != panel {n_full:,} — aborting")
        af_raw = t.alt_freq.to_numpy(dtype=np.float64)
        nc_raw = t.n_called.to_numpy()
        af = encode_af(af_raw)                              # u16, 4-decimal
        nc = np.clip(nc_raw, 0, 255).astype(np.uint8)
        p = _paths(d, name)
        # write to .tmp then rename so a killed task never leaves a partial vector
        for path, arr in ((p["af_nonsnp"], af[nonsnp]), (p["af_snp"], af[snp]),
                          (p["nc_nonsnp"], nc[nonsnp]), (p["nc_snp"], nc[snp])):
            np.save(path + ".tmp.npy", arr)
            os.replace(path + ".tmp.npy", path)
        print(f"  [{j}] {name} -> nonsnp {n_nonsnp:,} / snp {n_snp:,}", flush=True)
    print(f"convert [{s},{end}) complete")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("init", "convert"):
        p = sub.add_parser(name)
        p.add_argument("--out-dir", required=True)
        p.add_argument("--base", default=OUT)
        if name == "convert":
            p.add_argument("--start", type=int, required=True)
            p.add_argument("--count", type=int, required=True)
    args = ap.parse_args()
    {"init": cmd_init, "convert": cmd_convert}[args.cmd](args)


if __name__ == "__main__":
    main()
