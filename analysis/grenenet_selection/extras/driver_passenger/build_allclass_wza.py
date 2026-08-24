#!/usr/bin/env python
"""All-class WZA on clq0.9 LD blocks — step 1 of the driver-vs-passenger plan.

Takes the per-record kendall GEA outputs for the three classes (snp / sv /
smallindel), RE-STAMPS each record with its BigLD clq0.9 block (lib.assign_clq_blocks,
finer & tighter-LD than the phase-1 hapFIRE map that build_class_matrices used),
POOLS all three classes into ONE per-record table, and runs the canonical Booker
WZA per clq0.9 block.

Design decisions (this session):
  * all-class pooled screen (not per-class): the block is nominated on TOTAL,
    class-agnostic evidence -> the within-block conditional test then adjudicates
    which class drives. Screening on SV-only WZA would be circular.
  * clq0.9 blocks: tight LD islands (median ~223 bp / 7 var) -> ~single haplotype
    -> a clean unit for the later "does the SV beat the best SNP" fine-map. Gap
    records (not in any island) are dropped.
  * deg-2 WZA (canonical): clq0.9 removes the huge (9k-SNP) windows that made
    phase-1's deg-2 SNP-number-correction NaN; we VERIFY that here (NaN / neg-SD /
    fit-RMSE deg2 vs deg7) rather than defaulting to the deg-7 workaround.

Outputs (--out, default analysis/grenenet_selection/extras/driver_passenger/results):
  allclass_{model}_gen{g}_{climate}.records.csv   pooled per-record table
  wza_allclass_{model}_gen{g}_{climate}_deg{2,7}.csv
  block_composition_{model}_gen{g}_{climate}.csv  per-block class counts + WZA Z_pVal

Usage:
  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY build_allclass_wza.py --model kendall --gen 9 --climate bio1
"""
from __future__ import annotations
import argparse, os, subprocess, sys, tempfile
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GEA_DIR = os.path.dirname(HERE)
sys.path.insert(0, GEA_DIR)
import lib

WZA = os.path.join(GEA_DIR, "wza_script.py")
PY = sys.executable
SV_MIN_BP = 50
CLASSES = ("snp", "smallindel", "sv")


def build_records(model: str, gen: int, climate: str, indir: str) -> pd.DataFrame:
    """Load the 3 per-class kendall outputs, re-stamp clq0.9 block, pool."""
    parts = []
    for cls in CLASSES:
        f = f"{indir}/{model}_{cls}_gen{gen}_{climate}.csv"
        if not os.path.exists(f):
            sys.exit(f"missing per-class model result: {f}")
        d = pd.read_csv(f)
        d["cls"] = cls
        parts.append(d)
        print(f"  loaded {cls:11s} {len(d):>9,} records", flush=True)
    df = pd.concat(parts, ignore_index=True)

    # RE-STAMP: overwrite the phase-1 hapFIRE `block` with the clq0.9 island id.
    clq = lib.assign_clq_blocks(np.asarray(df["chrom"], dtype=str),
                                np.asarray(df["pos"]), r2=0.9)
    df["block"] = clq
    n_gap = (df["block"].astype(str) == "").sum()
    df = df[df["block"].astype(str) != ""].copy()
    print(f"  clq0.9 re-stamp: dropped {n_gap:,} inter-island gap records; "
          f"{len(df):,} kept in {df['block'].nunique():,} blocks", flush=True)
    return df


def run_wza(in_csv: str, out_csv: str, poly_deg: int, min_snps: int = 2):
    cmd = [PY, WZA, "--correlations", in_csv, "--summary_stat", "pval",
           "--window", "block", "--MAF", "MAF", "--maf_filter", "0.05",
           "--min_snps", str(min_snps), "--sep", ",", "--retain", "chrom", "pos",
           "--poly_deg", str(poly_deg), "--min_entries", "40", "--roller", "50",
           "--output", out_csv]
    print("  $ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=os.path.dirname(out_csv) or ".")


def _rolling_support(wza_df: pd.DataFrame, roller=50, min_entries=40):
    """Reproduce adjust_WZA_with_spline's rolling SD-vs-meanSNP support curve."""
    s = wza_df[~wza_df["Z"].isnull()].sort_values("SNPs")
    var = s["Z"].rolling(window=roller, min_periods=min_entries).var()
    mask = ~var.isnull()
    sd = np.sqrt(var)[mask].to_numpy()
    mean_snp = s["SNPs"].rolling(window=roller, min_periods=min_entries).mean()[mask].to_numpy()
    return mean_snp, sd


def poly_diag(wza_df: pd.DataFrame, snp_grid: np.ndarray):
    """RMSE of the SD-vs-nSNP fit for deg2 & deg7 + how many windows would get a
    non-positive predicted SD (the pre-floor NaN cause) at each degree."""
    mean_snp, sd = _rolling_support(wza_df)
    out = {}
    for deg in (2, 7):
        if len(mean_snp) <= deg:
            out[deg] = dict(rmse=np.nan, neg_sd_windows=np.nan, n_support=len(mean_snp))
            continue
        coef = np.polyfit(mean_snp, sd, deg=deg)
        model = np.poly1d(coef)
        rmse = float(np.sqrt(np.mean((model(mean_snp) - sd) ** 2)))
        pred = model(snp_grid)
        out[deg] = dict(rmse=rmse, neg_sd_windows=int((pred <= 0).sum()),
                        n_support=len(mean_snp))
    return out


def block_composition(recs: pd.DataFrame, wza: pd.DataFrame) -> pd.DataFrame:
    """Per-block class counts joined to the WZA block-level Z_pVal."""
    g = recs.groupby("block")
    comp = pd.DataFrame({
        "n_rec": g.size(),
        "n_snp": g.apply(lambda x: (x["cls"] == "snp").sum()),
        "n_indel": g.apply(lambda x: (x["cls"] == "smallindel").sum()),
        "n_sv": g.apply(lambda x: (x["cls"] == "sv").sum()),
        "chrom": g["chrom"].first(),
        "start": g["pos"].min(),
        "end": g["pos"].max(),
    }).reset_index()
    comp["has_sv"] = comp["n_sv"] > 0
    comp["sv_excl"] = (comp["n_snp"] == 0) & (comp["n_indel"] == 0) & (comp["n_sv"] > 0)
    # wza_script emits the window-id column as "index" (retain path) or "gene"
    idcol = "index" if "index" in wza.columns else "gene"
    if "Z_pVal" in wza.columns:
        comp = comp.merge(wza[[idcol, "Z_pVal", "Z", "SNPs"]]
                          .rename(columns={idcol: "block"}), on="block", how="left")
    return comp


def report_composition(comp: pd.DataFrame):
    tot = len(comp)
    scored = comp["Z_pVal"].notna() if "Z_pVal" in comp else pd.Series(False, index=comp.index)
    print("\n" + "=" * 68)
    print("BLOCK COMPOSITION (clq0.9, all-class pooled)")
    print("=" * 68)
    print(f"  blocks with >=1 record        : {tot:,}")
    print(f"  blocks WZA-scored (>=2 rec)   : {int(scored.sum()):,}")
    print(f"  blocks with >=1 SV            : {int(comp['has_sv'].sum()):,} "
          f"({100*comp['has_sv'].mean():.1f}%)")
    sv_excl = comp["sv_excl"]
    print(f"  SV-exclusive blocks           : {int(sv_excl.sum()):,} "
          f"({100*sv_excl.mean():.1f}%)")
    print(f"    of which singleton (1 SV)   : {int((sv_excl & (comp['n_rec']==1)).sum()):,}"
          "   -> route to window tag-SNP test, not WZA")
    print(f"    SV-exclusive & WZA-scored   : {int((sv_excl & scored).sum()):,}")
    print("  per-block size (n_rec): "
          + ", ".join(f"{q}={comp['n_rec'].quantile(p):.0f}"
                      for q, p in [("p50", .5), ("p90", .9), ("p99", .99), ("max", 1.0)]))
    if scored.any():
        s = comp[scored].copy()
        for k in (50, 200, 1000):
            top = s.nsmallest(k, "Z_pVal")
            print(f"  top {k:>4} WZA blocks: {int(top['has_sv'].sum()):>4} contain an SV "
                  f"({100*top['has_sv'].mean():.0f}%); "
                  f"{int((top['n_snp']==0).sum())} SNP-free")
    print("=" * 68)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="kendall")
    ap.add_argument("--gen", type=int, default=9)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--indir", default=None,
                    help="dir of per-class model CSVs (default ../phase1_replication/<model>)")
    ap.add_argument("--out", default=f"{lib.GEA}/extras/driver_passenger/results")
    ap.add_argument("--min-snps", type=int, default=2)
    args = ap.parse_args()

    indir = args.indir or f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/{args.model}"
    os.makedirs(args.out, exist_ok=True)
    stem = f"allclass_{args.model}_gen{args.gen}_{args.climate}"

    print(f"[1/3] pooling per-class records ({args.model} gen{args.gen} {args.climate})", flush=True)
    recs = build_records(args.model, args.gen, args.climate, indir)
    rec_csv = f"{args.out}/{stem}.records.csv"
    recs.to_csv(rec_csv, index=False)
    print(f"  -> {rec_csv}", flush=True)

    # clean input for WZA (drop any NaN pval / empty block already handled)
    clean = recs[recs["pval"].notna()].copy()
    with tempfile.NamedTemporaryFile("w", suffix=".csv", dir=args.out, delete=False) as tf:
        clean.to_csv(tf.name, index=False)
        tmp = tf.name

    wza_out = {}
    try:
        for deg in (2, 7):
            o = f"{args.out}/wza_{stem}_deg{deg}.csv"
            print(f"\n[2/3] WZA deg-{deg}", flush=True)
            run_wza(tmp, o, poly_deg=deg, min_snps=args.min_snps)
            wza_out[deg] = pd.read_csv(o)
    finally:
        os.unlink(tmp)

    # deg-2 verification: NaN Z, neg-SD-pre-floor, fit RMSE
    w2 = wza_out[2]
    snp_grid = w2["SNPs"].to_numpy()
    diag = poly_diag(w2, snp_grid)
    print("\n" + "=" * 68)
    print("WZA SNP-NUMBER-CORRECTION DIAGNOSTIC (deg-2 verification)")
    print("=" * 68)
    print(f"  windows scored              : {len(w2):,}")
    print(f"  NaN Z (raw weiZ)            : {int(w2['Z'].isna().sum()):,}")
    print(f"  NaN Z_pVal (post-floor)     : {int(w2['Z_pVal'].isna().sum()):,}")
    print(f"  max SNPs/window             : {int(w2['SNPs'].max()):,} "
          f"(windows >2000: {int((w2['SNPs']>2000).sum())})")
    print(f"  support-curve points        : {diag[2]['n_support']:,}")
    for deg in (2, 7):
        d = diag[deg]
        print(f"  deg-{deg}: SD-fit RMSE={d['rmse']:.3f}  "
              f"neg/zero-SD windows (pre-floor)={d['neg_sd_windows']}")
    print("=" * 68)

    print("\n[3/3] block composition", flush=True)
    comp = block_composition(recs, wza_out[2])
    comp_csv = f"{args.out}/block_composition_{args.model}_gen{args.gen}_{args.climate}.csv"
    comp.to_csv(comp_csv, index=False)
    report_composition(comp)
    print(f"\n  -> {comp_csv}", flush=True)


if __name__ == "__main__":
    main()
