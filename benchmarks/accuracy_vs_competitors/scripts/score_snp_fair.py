#!/usr/bin/env python3
"""Allele-aware (fair) SNP-AF scoring of kMate / hapFIRE / vg vs sim truth.

WHY THIS EXISTS
---------------
The earlier scorers (benchmarks/score_competitor_snp.py and
benchmarks/accuracy_vs_competitors/scripts/score_competitors.py) join every tool
to the truth on (chrom,pos) ONLY. The sim `recomb_truth.tsv.gz` and the kMate
`_global.tsv` carry NO ref/alt bases -- only ref_len/alt_len -- and the truth
has ~45k positions with MORE THAN ONE biallelic-SNP record (multiallelic sites).
A position-only join therefore pairs the WRONG allele between a tool's estimate
and the truth at those positions, deflating R^2 (e.g. single-founder pool gave a
spurious kMate R^2=0.937 instead of the physically-required ~1.0).

THE FIX
-------
Score every tool on the SHARED SNP PANEL -- the bcftools-norm'd, one-SNP-per-pos,
biallelic panel (shared_snps_p80_Chr1.vcf.gz) with known REF/ALT. Empirically,
ALL of the truth's multiallelic/duplicate positions were dropped when the panel
was built, so AFTER restricting truth (and kMate, and the competitors) to the
panel positions, (chrom,pos) is a UNIQUE, unambiguous key for every tool -- one
SNP per position, one truth row per position, one kMate row per position. The
position-only join that hapFIRE/vg require (they emit only chrom,pos,freq) is
then correct by construction, and kMate's allele is pinned to the same panel SNP.

This was verified: on cov10_n1_g0_s42 (single founder) the panel join yields
kMate R^2 ~= 1.0, vs 0.937 from the buggy join.

SCORING BASIS
-------------
truth_af is the MAR convention (AF divided by *called* founders). At well-informed
sites (info>=0.9) MAR == physical AF, so info>=0.9 is the fair, convention-free
basis for comparing a model-based caller (kMate) against read-based tools
(hapFIRE/vg). At low info, scoring vs truth_af FAVORS kMate (closed loop) -- see
the stratification mode and the caveat in --help.
"""
import argparse, gzip, sys
from pathlib import Path
import numpy as np
import pandas as pd


# --------------------------------------------------------------------------- #
def norm_chrom(s):
    return (s.astype(str).str.replace("Chr", "", regex=False)
                          .str.replace("chr", "", regex=False))


def load_panel_positions(vcf_gz):
    """Return set of int positions for the shared (one-SNP-per-pos) panel."""
    pos = set()
    with gzip.open(vcf_gz, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t", 3)
            pos.add(int(c[1]))
    return pos


def load_truth(path, panel_pos):
    t = pd.read_csv(path, sep="\t")
    t = t[(t.ref_len == 1) & (t.alt_len == 1)].copy()          # SNP truth rows
    t = t[t["pos"].isin(panel_pos)].copy()                     # -> panel SNP set
    # after panel restriction there is exactly one truth SNP row per position
    dup = t["pos"].duplicated().sum()
    if dup:
        sys.stderr.write(f"WARNING: {dup} duplicate truth positions survived "
                         f"panel restriction; keeping first\n")
        t = t.drop_duplicates("pos", keep="first")
    return t


def load_kmate(path, panel_pos):
    d = pd.read_csv(path, sep="\t")
    d = d[(d.ref_len == 1) & (d.alt_len == 1)].copy()          # SNP rows
    d = d[d["pos"].isin(panel_pos)].copy()                     # -> panel SNP set
    d = d.drop_duplicates("pos", keep="first")
    return d[["pos", "alt_freq"]].rename(columns={"alt_freq": "est"})


def load_3col(path, panel_pos):
    """hapFIRE / vg: chrom, pos, freq (no header)."""
    d = pd.read_csv(path, sep="\t", header=None, names=["chrom", "pos", "est"],
                    na_values=["nan"])
    d["pos"] = pd.to_numeric(d["pos"], errors="coerce").astype("Int64")
    d["est"] = pd.to_numeric(d["est"], errors="coerce")
    d = d[d["pos"].isin(panel_pos)].copy()
    d = d.drop_duplicates("pos", keep="first")
    return d[["pos", "est"]]


def metrics(truth, est):
    truth = np.asarray(truth, float)
    est = np.asarray(est, float)
    m = np.isfinite(truth) & np.isfinite(est)
    truth, est = truth[m], est[m]
    if len(truth) == 0:
        return dict(n=0, MAE=np.nan, RMSE=np.nan, R2=np.nan, pearson_r=np.nan)
    err = est - truth
    mae = float(np.abs(err).mean())
    rmse = float(np.sqrt((err ** 2).mean()))
    ss_res = float((err ** 2).sum())
    ss_tot = float(((truth - truth.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    pear = float(np.corrcoef(truth, est)[0, 1]) if len(truth) > 1 else np.nan
    return dict(n=int(len(truth)), MAE=mae, RMSE=rmse, R2=r2, pearson_r=pear)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", required=True)
    ap.add_argument("--panel", required=True, help="shared-SNP-panel VCF(.gz)")
    ap.add_argument("--truth", required=True, help="recomb_truth.tsv.gz")
    ap.add_argument("--kmate", help="kMate <pool>_global.tsv")
    ap.add_argument("--hapfire", help="hapFIRE _snp_frequency.txt (chrom pos freq)")
    ap.add_argument("--vg", help="vg _snp_frequency.txt (chrom pos freq)")
    ap.add_argument("--info-min", type=float, default=0.9,
                    help="score only truth SNPs with info>=this (fair basis where "
                         "MAR==physical). NOTE truth_af is MAR; at low info scoring "
                         "vs it FAVORS kMate (closed loop). default 0.9")
    ap.add_argument("--out", help="append tidy rows (pool,tool,n,MAE,RMSE,R2,pearson_r)")
    ap.add_argument("--stratify", action="store_true",
                    help="instead of a single cut, report metrics per info bin "
                         "[0,.5) [.5,.9) [.9,.99) [.99,1] for each tool")
    a = ap.parse_args()

    panel_pos = load_panel_positions(a.panel)
    tr = load_truth(a.truth, panel_pos)
    sys.stderr.write(f"[{a.pool}] panel SNPs={len(panel_pos):,}  "
                     f"truth-SNPs-on-panel={len(tr):,}\n")

    tools = {}
    if a.kmate:
        tools["kMate"] = load_kmate(a.kmate, panel_pos)
    if a.hapfire:
        tools["hapFIRE"] = load_3col(a.hapfire, panel_pos)
    if a.vg:
        tools["vg-giraffe"] = load_3col(a.vg, panel_pos)

    rows = []
    if a.stratify:
        bins = [(-0.01, 0.5), (0.5, 0.9), (0.9, 0.99), (0.99, 1.0001)]
        for name, est in tools.items():
            j = tr.merge(est, on="pos", how="inner")
            for lo, hi in bins:
                sub = j[(j["info"] >= lo) & (j["info"] < hi)]
                mt = metrics(sub["truth_af"].values, sub["est"].values)
                lab = f"[{max(lo,0):.2f},{hi if hi<=1 else 1.0:.2f})"
                rows.append({"pool": a.pool, "tool": name, "info_bin": lab, **mt})
                sys.stderr.write(f"  {name:11s} info {lab}: n={mt['n']:>8,} "
                                 f"R2={mt['R2'] if mt['R2']==mt['R2'] else float('nan'):.4f} "
                                 f"MAE={mt['MAE']:.4f}\n")
        out = pd.DataFrame(rows)[["pool", "tool", "info_bin", "n", "MAE",
                                  "RMSE", "R2", "pearson_r"]]
        print(out.to_string(index=False))
        if a.out:
            tp = Path(a.out)
            out.to_csv(tp, sep="\t", mode="a", header=not tp.exists(), index=False)
        return

    trc = tr[tr["info"] >= a.info_min] if a.info_min > 0 else tr
    sys.stderr.write(f"  info>={a.info_min}: {len(trc):,}/{len(tr):,} truth SNPs\n")
    for name, est in tools.items():
        j = trc.merge(est, on="pos", how="inner")
        mt = metrics(j["truth_af"].values, j["est"].values)
        rows.append({"pool": a.pool, "tool": name, **mt})
        sys.stderr.write(f"  {name:11s}: n={mt['n']:>8,}  MAE={mt['MAE']:.4f}  "
                         f"RMSE={mt['RMSE']:.4f}  R2={mt['R2']:.4f}  "
                         f"r={mt['pearson_r']:.4f}\n")
    out = pd.DataFrame(rows)[["pool", "tool", "n", "MAE", "RMSE", "R2", "pearson_r"]]
    print(out.to_string(index=False))
    if a.out:
        tp = Path(a.out)
        out.to_csv(tp, sep="\t", mode="a", header=not tp.exists(), index=False)


if __name__ == "__main__":
    main()
