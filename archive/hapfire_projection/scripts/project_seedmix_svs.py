#!/usr/bin/env python3
"""
Project hapFIRE seedmix ecotype frequencies onto SV ALT frequencies via the
founder × SV genotype matrix; compare to freqk seedmix per-SV outputs.

Inputs:
  --hapfire-dir   directory of hapFIRE seed_mix outputs (s{1..8}_ecotype_frequency.txt)
  --founder-mat   data/founder_sv_matrix.parquet (rows=SV records, cols=ecotypes)
  --founder-meta  data/founder_sv_meta.parquet  (chrom, pos, alt_idx, var_type, sv_size, ac_in_panel)
  --freqk-svs     freqk seedmix_p_svs.parquet (chrom, pos, SEEDMIX_S1..S8)
                  OR seedmix_p_wide.parquet (with alt_idx column)
  --variants      variants.tsv (chrom, pos, alt_idx, var_type, sv_size) — used to align indexes

Outputs:
  hapfire_proj_svs_seedmix.parquet  (chrom, pos, alt_idx, hapfire_S1..S8)
  hapfire_vs_freqk_svs_seedmix_AF.csv.gz  (per-row tall comparison)
  hapfire_vs_freqk_svs_seedmix_per_sample_corr.tsv (per-sample r/r2)
"""
import argparse
import gzip
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

SAMPLE_IDS = [f"SEEDMIX_S{i}" for i in range(1, 9)]


def load_hapfire_ecotype(path: Path) -> dict:
    df = pd.read_csv(path, sep="\t", header=None, names=["ecotype", "freq"], dtype={"ecotype": str})
    return dict(zip(df["ecotype"], df["freq"].astype("float64")))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hapfire-dir",
                    default="/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_frequencies/seed_mix")
    ap.add_argument("--founder-mat", required=True, type=Path)
    ap.add_argument("--founder-meta", required=True, type=Path)
    ap.add_argument("--freqk-wide",
                    default="/home/tbellagio/scratch/freqk_gr/results/seedmix_p_wide.parquet")
    ap.add_argument("--variants",
                    default="/home/tbellagio/scratch/freqk_gr/results/variants.tsv")
    ap.add_argument("--out-proj", required=True, type=Path)
    ap.add_argument("--out-tall", required=True, type=Path)
    ap.add_argument("--out-corr", required=True, type=Path)
    args = ap.parse_args()

    t0 = time.time()
    print(f"[+{time.time()-t0:.0f}s] Loading founder × SV matrix from {args.founder_mat}")
    G = pd.read_parquet(args.founder_mat)
    meta = pd.read_parquet(args.founder_meta)
    print(f"  G shape: {G.shape}; meta rows: {len(meta)}")

    ecotype_cols = [c for c in G.columns if c not in ("chrom", "pos", "alt_idx")]
    G_keys = G[["chrom", "pos", "alt_idx"]].copy()
    G_keys["chrom"] = G_keys["chrom"].astype(str)
    G_keys["pos"] = G_keys["pos"].astype(np.int32)
    G_keys["alt_idx"] = G_keys["alt_idx"].astype(np.int8)
    G_mat = G[ecotype_cols].to_numpy(dtype=np.float32)
    print(f"  {len(ecotype_cols)} ecotype columns; {G_mat.shape[0]} SV records")
    print(f"  ecotypes head: {ecotype_cols[:5]}")

    # collapse duplicate ecotype columns (e.g. multiple Assembly_IDs → same Accession_ID).
    # If two Assembly_IDs both report 1 for an ALT, the column is still 1 (not a sum).
    df_e = pd.DataFrame(G_mat, columns=ecotype_cols)
    G_collapsed = df_e.groupby(df_e.columns, axis=1).max()
    G_mat = G_collapsed.to_numpy(dtype=np.float32)
    ecotype_cols = G_collapsed.columns.tolist()
    print(f"  after dedup of duplicate ecotype cols: {len(ecotype_cols)}")

    print(f"\n[+{time.time()-t0:.0f}s] Loading hapFIRE seedmix ecotype freqs")
    hapfire = {}
    for sid in SAMPLE_IDS:
        n = sid.split("_S")[1]
        path = Path(args.hapfire_dir) / f"s{n}_ecotype_frequency.txt"
        hapfire[sid] = load_hapfire_ecotype(path)
        print(f"  {sid}: {len(hapfire[sid])} ecotypes; sum={sum(hapfire[sid].values()):.4f}")

    # Build f_ecotype matrix (8 samples × len(ecotype_cols))
    F = np.zeros((len(SAMPLE_IDS), len(ecotype_cols)), dtype=np.float32)
    missing_per_sample = {}
    for i, sid in enumerate(SAMPLE_IDS):
        h = hapfire[sid]
        missing = []
        for j, e in enumerate(ecotype_cols):
            if e in h:
                F[i, j] = h[e]
            else:
                missing.append(e)
        missing_per_sample[sid] = missing
    if any(missing_per_sample.values()):
        sample0_miss = missing_per_sample[SAMPLE_IDS[0]]
        print(f"  missing in hapFIRE panel: {len(sample0_miss)} ecotypes (e.g. {sample0_miss[:5]})")
    print(f"  F shape: {F.shape}; row sums (mass on mapped 80): {F.sum(axis=1).round(3)}")

    # Project: f_proj_per_SV[s, v] = F[s, :] @ G_mat[v, :]  ⇔  P = F @ G_mat.T
    print(f"\n[+{time.time()-t0:.0f}s] Projecting: F @ G.T")
    P = F @ G_mat.T  # (8, n_sv)
    print(f"  P shape: {P.shape}  (samples × SV records)")

    proj = G_keys.copy()
    for i, sid in enumerate(SAMPLE_IDS):
        proj[sid] = P[i].astype("float32")
    proj["chrom"] = proj["chrom"].astype("category")
    proj.to_parquet(args.out_proj, compression="zstd")
    print(f"  wrote: {args.out_proj}  ({args.out_proj.stat().st_size/1e6:.1f} MB)")

    # Per-sample mass on mapped 80 founders (call this 'panel_coverage')
    panel_cov = F.sum(axis=1)
    print("\nPer-sample hapFIRE mass on mapped 80 founders:")
    for sid, m in zip(SAMPLE_IDS, panel_cov):
        print(f"  {sid}: {m:.4f}  (=> SV freq projections are scaled by this)")

    # Load freqk seedmix wide and join on (chrom, pos, alt_idx)
    print(f"\n[+{time.time()-t0:.0f}s] Loading freqk seedmix wide: {args.freqk_wide}")
    freqk = pd.read_parquet(args.freqk_wide)
    freqk["chrom"] = freqk["chrom"].astype(str)
    freqk["pos"] = freqk["pos"].astype(np.int32)
    freqk["alt_idx"] = freqk["alt_idx"].astype(np.int8)
    print(f"  freqk wide rows: {len(freqk):,}")

    # restrict freqk to SV rows (matching our G_keys)
    proj["chrom"] = proj["chrom"].astype(str)
    merged = freqk.merge(proj, on=["chrom", "pos", "alt_idx"], suffixes=("_freqk", "_hapf"))
    print(f"  merged on (chrom, pos, alt_idx): {len(merged):,} SVs in both")

    # also tag with var_type via meta
    meta["chrom"] = meta["chrom"].astype(str)
    meta_keys = meta[["chrom", "pos", "alt_idx", "var_type", "sv_size", "ac_in_panel"]].copy()
    merged = merged.merge(meta_keys, on=["chrom", "pos", "alt_idx"])
    print(f"  with metadata: {len(merged):,}")

    # Per-sample correlation
    corr_rows = []
    for sid in SAMPLE_IDS:
        x = merged[f"{sid}_freqk"].to_numpy(dtype="float64")
        y = merged[f"{sid}_hapf"].to_numpy(dtype="float64")
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() > 1:
            r = float(np.corrcoef(x[ok], y[ok])[0, 1])
            mae = float(np.abs(x[ok] - y[ok]).mean())
            rmse = float(np.sqrt(((x[ok] - y[ok])**2).mean()))
        else:
            r = mae = rmse = float("nan")
        corr_rows.append({"sample_id": sid, "n_sv": int(ok.sum()),
                          "pearson_r": r, "r_squared": r*r,
                          "mae": mae, "rmse": rmse,
                          "panel_coverage": float(panel_cov[SAMPLE_IDS.index(sid)])})
    corr = pd.DataFrame(corr_rows)
    corr.to_csv(args.out_corr, sep="\t", index=False)
    print(f"\n  wrote: {args.out_corr}")
    print(corr.round(4).to_string(index=False))

    # write tall comparison
    tall_rows = []
    for sid in SAMPLE_IDS:
        sub = merged[["chrom", "pos", "alt_idx", "var_type", "sv_size", "ac_in_panel",
                      f"{sid}_freqk", f"{sid}_hapf"]].rename(
            columns={f"{sid}_freqk": "freqk_af", f"{sid}_hapf": "hapfire_proj_af"}
        )
        sub.insert(0, "sample_id", sid)
        tall_rows.append(sub)
    tall = pd.concat(tall_rows, ignore_index=True)
    tall.to_csv(args.out_tall, sep="\t", index=False, compression="gzip", float_format="%.4f")
    print(f"  wrote: {args.out_tall}  ({args.out_tall.stat().st_size/1e6:.1f} MB)")

    # By var_type breakdown for the first sample
    print("\nPer var_type, sample S1 stats (n / mean abs error / r):")
    s = "SEEDMIX_S1"
    for vt in sorted(merged["var_type"].unique()):
        sub = merged[merged["var_type"] == vt]
        x = sub[f"{s}_freqk"].to_numpy(dtype="float64")
        y = sub[f"{s}_hapf"].to_numpy(dtype="float64")
        ok = np.isfinite(x) & np.isfinite(y)
        n = int(ok.sum())
        if n > 1:
            r = float(np.corrcoef(x[ok], y[ok])[0, 1])
            mae = float(np.abs(x[ok] - y[ok]).mean())
            print(f"  {vt}: n={n:7d}  mae={mae:.4f}  r={r:.4f}")
        else:
            print(f"  {vt}: n={n}")


if __name__ == "__main__":
    main()
