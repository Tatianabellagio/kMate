#!/usr/bin/env python3
"""Score kMate (global+block), hapFIRE, vg(-SNP, -SV) into ONE 4-tool benchmark
table for the paper RMSE figure. Emits a `basis` column with BOTH conventions:

  basis=allrec     -> all records, MAR truth_af (the existing kMate-table basis;
                      the design doc warns this flatters kMate / sinks read-based tools)
  basis=fullcalled -> convention-free subset where MAR == physical:
                        SNP: info >= --fullcalled-info (0.99)
                        SV : panel n_called == F (all founders genotyped)

Reuses the VALIDATED joins so SV R^2 isn't crushed by a bad key:
  SNP -> score_snp_fair (shared-SNP-panel (chrom,pos) join vs recomb_truth)
  SV  -> svidx (row index over the SV subset of recomb_truth, in var_pa-meta order),
         the same key score_sv.py uses. vg-SV est is already (svidx, est).

recomb_truth carries valid SV truth_af for ALL generations (g0/g1/g3 mosaic), so
this path works across the whole grid (build_truth_sv.py is g0-only -- not used).

One invocation per pool; appends rows to --table.
Cols: tool mode panel coverage n_founders generation mating selection seed
      var_class basis n MAE RMSE R2 pearson_r est_file
"""
import argparse, re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import score_snp_fair as ssf   # load_panel_positions, load_truth, load_kmate, load_3col, metrics

COLS = ["tool", "mode", "panel", "coverage", "n_founders", "generation", "mating",
        "selection", "seed", "var_class", "basis", "n", "MAE", "RMSE", "R2",
        "pearson_r", "est_file"]


def parse_pool(pool):
    def g(pat, cast=str):
        m = re.search(pat, pool)
        return cast(m.group(1)) if m else None
    sel = "dom500nr" if "dom500nr" in pool else ("dom500" if "dom500" in pool else "balanced")
    return dict(panel="p80" if "p80" in pool else ("p231" if "p231" in pool else "NA"),
                coverage=g(r"cov(\d+)", int), n_founders=g(r"_n(\d+)", int),
                generation=g(r"_g(\d+)", int),
                mating="self97" if "self97" in pool else "outcross",
                selection=sel, seed=g(r"_s(\d+)", int))


def emit(rows, meta, tool, mode, vclass, basis, mt, est_file):
    rows.append({**meta, "tool": tool, "mode": mode, "var_class": vclass,
                 "basis": basis, **mt, "est_file": est_file})
    r2 = mt["R2"]
    sys.stderr.write(f"  {tool:11s}/{mode:6s} {vclass:3s} {basis:10s}: "
                     f"n={mt['n']:>8,} RMSE={mt['RMSE'] if mt['RMSE']==mt['RMSE'] else float('nan'):.4f} "
                     f"R2={r2 if r2==r2 else float('nan'):.4f}\n")


def score_snp(rows, meta, args, panel_pos):
    tr = ssf.load_truth(args.truth, panel_pos)          # SNP truth rows on shared panel
    ests = {}   # (tool, mode) -> (est_df, est_path)
    if args.kmate_global:
        ests[("kMate", "global")] = (ssf.load_kmate(args.kmate_global, panel_pos), args.kmate_global)
    if args.kmate_block:
        ests[("kMate", "block")] = (ssf.load_kmate(args.kmate_block, panel_pos), args.kmate_block)
    if args.hapfire:
        ests[("hapFIRE", "default")] = (ssf.load_3col(args.hapfire, panel_pos), args.hapfire)
    if args.vg_snp:
        ests[("vg", "default")] = (ssf.load_3col(args.vg_snp, panel_pos), args.vg_snp)
    for basis, trc in [("allrec", tr), ("fullcalled", tr[tr["info"] >= args.fullcalled_info])]:
        for (tool, mode), (est, path) in ests.items():
            j = trc.merge(est, on="pos", how="inner")
            mt = ssf.metrics(j["truth_af"].values, j["est"].values)
            emit(rows, meta, tool, mode, "SNP", basis, mt, Path(path).name)


def score_sv(rows, meta, args):
    tr = pd.read_csv(args.truth, sep="\t")
    m = np.load(args.var_meta, allow_pickle=True)
    nrec_meta = len(m["pos"])
    if len(tr) != nrec_meta:
        sys.exit(f"[sv] FATAL: recomb_truth rows {len(tr):,} != var_pa meta {nrec_meta:,} "
                 "-- svidx alignment unsafe.")
    svmask = (np.abs(m["alt_len"].astype(int) - m["ref_len"].astype(int)) >= args.svlen)
    vc = sp.load_npz(args.var_called).tocsr()
    n_called_sv = np.asarray(vc[:, svmask].sum(axis=0)).ravel().astype(int)   # panel n_called per SV
    F = vc.shape[0]
    trV = tr[svmask].reset_index(drop=True).copy()
    trV["svidx"] = np.arange(len(trV))
    trV["n_called"] = n_called_sv
    assert len(trV) == len(n_called_sv)

    def load_kmate_sv(p):
        d = pd.read_csv(p, sep="\t")
        d = d[(np.abs(d.alt_len - d.ref_len) >= args.svlen)].reset_index(drop=True)
        d["svidx"] = np.arange(len(d))
        return d[["svidx", "alt_freq"]].rename(columns={"alt_freq": "est"})

    ests = {}
    if args.kmate_global:
        ests[("kMate", "global")] = load_kmate_sv(args.kmate_global)
    if args.kmate_block:
        ests[("kMate", "block")] = load_kmate_sv(args.kmate_block)
    if args.vg_sv:
        ests[("vg", "default")] = pd.read_csv(args.vg_sv, sep="\t")[["svidx", "est"]]

    for basis, sub in [("allrec", trV), ("fullcalled", trV[trV["n_called"] == F])]:
        for (tool, mode), est in ests.items():
            j = sub.merge(est, on="svidx", how="inner")
            mt = ssf.metrics(j["truth_af"].values, j["est"].values)
            ef = Path(args.vg_sv).name if tool == "vg" else Path(getattr(args, "kmate_" + mode)).name
            emit(rows, meta, tool, mode, "SV", basis, mt, ef)
        # hapFIRE: explicit blank SV row (no native SV estimation, design doc §3.1)
        emit(rows, meta, "hapFIRE", "default", "SV", basis,
             dict(n=0, MAE=np.nan, RMSE=np.nan, R2=np.nan, pearson_r=np.nan), "NA")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", required=True)
    ap.add_argument("--truth", required=True, help="recomb_truth.tsv.gz")
    ap.add_argument("--snp-panel", required=True, help="shared SNP panel VCF.gz")
    ap.add_argument("--var-meta", required=True)
    ap.add_argument("--var-called", required=True)
    ap.add_argument("--kmate-global"); ap.add_argument("--kmate-block")
    ap.add_argument("--hapfire"); ap.add_argument("--vg-snp"); ap.add_argument("--vg-sv")
    ap.add_argument("--svlen", type=int, default=50)
    ap.add_argument("--fullcalled-info", type=float, default=0.99)
    ap.add_argument("--table", required=True)
    args = ap.parse_args()

    meta = parse_pool(args.pool)
    sys.stderr.write(f"[{args.pool}] {meta}\n")
    rows = []
    panel_pos = ssf.load_panel_positions(args.snp_panel)
    score_snp(rows, meta, args, panel_pos)
    score_sv(rows, meta, args)

    out = pd.DataFrame(rows)[COLS]
    tp = Path(args.table)
    out.to_csv(tp, sep="\t", mode="a", header=not tp.exists(), index=False)
    sys.stderr.write(f"  -> appended {len(out)} rows to {tp}\n")


if __name__ == "__main__":
    main()
