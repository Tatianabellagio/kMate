#!/usr/bin/env python
"""What does a founder call-rate floor do to the GEA hit list?

Two things, in order of cost:

  --mode posthoc   apply the floor to the EXISTING per-record p-values and recompute the
                   per-(class x axis) Bonferroni threshold on the retained set. This is not
                   the same as re-running -- LFMM fits its 16 latent factors to whatever
                   records are in the matrix, so removing 12% of them shifts every p-value
                   slightly -- but it is a good forecast of the re-run, and it costs minutes
                   instead of 88 model fits.

  --mode pool      the same floor applied to the candidate pool
                   (genes/convergence/results/evidence_matrix.csv), so we can see which
                   candidates and which lines of evidence survive.

Why 0.50: within matched MAF strata, records called in >=50% of founders are calibrated like
fully-called ones (lambda 1.7-1.8), while the 78-only stratum runs lambda 2.1-2.3 and about
twice the hit rate -- see multiaxis/callrate_calibration.py. 78/231 = 0.338 falls below the
floor; 0.90 would additionally discard the well-behaved 0.50-0.90 band, which is most of the
SV layer.

env: kmate. Compute node.
"""
from __future__ import annotations
import argparse
import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..")))
import lib                                                        # noqa: E402

MULTI = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis"
WZAIN = f"{MULTI}/_superseded_pre_cr50/wza_in_clq09_tile"   # the PRE-floor run, on purpose: this script measures what the 0.50 call-rate floor removes
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1", "pc2", "pc3"]
FLOOR = 0.50


def call_rates(kind: str) -> pd.DataFrame:
    idx = np.load(f"{lib.AF_STORE}/index_{kind}.npz", allow_pickle=True)
    ch = idx["chrom"].astype(str)
    rl = idx["ref_len"].astype(np.int64)
    al = idx["alt_len"].astype(np.int64)
    n_cal = np.zeros(len(ch), np.int64)
    off = 0
    for c in ("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"):
        n = int((ch == c).sum())
        base = f"{lib.PROJ}/panel/arch3/{c.lower()}/var_pa_231_arch3_{c.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        prl = meta["ref_len"].astype(np.int64)
        pal = meta["alt_len"].astype(np.int64)
        issnp = (prl == 1) & (pal == 1)
        panel = np.where(issnp if kind == "snp" else ~issnp)[0]
        assert panel.size == n and np.array_equal(prl[panel], rl[off:off + n])
        vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        n_cal[off:off + n] = np.asarray((vc[:, panel] > 0).sum(0)).ravel()
        off += n
    return pd.DataFrame(dict(chrom=ch, pos=idx["pos"], ref_len=rl, alt_len=al,
                             call_rate=n_cal / 231.0))


def posthoc(classes: list[str], floor: float) -> None:
    rows = []
    for cls in classes:
        CR = call_rates("snp" if cls == "snp" else "nonsnp")
        for ax in AXES:
            f = f"{WZAIN}/lfmm_{cls}_gen9_{ax}.csv"
            if not os.path.exists(f):
                continue
            d = pd.read_csv(f, usecols=["chrom", "pos", "ref_len", "alt_len", "MAF", "pval"])
            d = d[d.MAF > 0.05]
            thr_old = 0.05 / len(d)
            old_hits = int((d.pval < thr_old).sum())
            d = d.merge(CR, on=["chrom", "pos", "ref_len", "alt_len"], how="left")
            k = d.call_rate >= floor
            sub = d[k]
            thr_new = 0.05 / len(sub)
            rows.append(dict(cls=cls, axis=ax, n_old=len(d), n_new=int(k.sum()),
                             hits_old=old_hits, hits_new=int((sub.pval < thr_new).sum()),
                             hits_kept_of_old=int(((d.pval < thr_old) & k).sum()),
                             bar_old=-np.log10(thr_old), bar_new=-np.log10(thr_new)))
            print(f"  {cls:11s} {ax:5s} {old_hits:>5,} -> {rows[-1]['hits_new']:>5,}", flush=True)
    D = pd.DataFrame(rows)
    D.to_csv(f"{MULTI}/callrate_floor_posthoc.csv", index=False)
    g = D.groupby("cls").agg(tests_old=("n_old", "sum"), tests_new=("n_new", "sum"),
                             hits_old=("hits_old", "sum"), hits_new=("hits_new", "sum"),
                             hits_of_old_kept=("hits_kept_of_old", "sum"))
    g["pct_tests_kept"] = (100 * g.tests_new / g.tests_old).round(1)
    g["pct_hits_kept"] = (100 * g.hits_new / g.hits_old).round(1)
    pd.set_option("display.width", 200)
    print(f"\nforecast of a call-rate >= {floor} floor (per class, summed over 22 axes)\n")
    print(g.to_string())
    print(f"\nwrote {MULTI}/callrate_floor_posthoc.csv")


def pool(floor: float) -> None:
    E = pd.read_csv(f"{lib.GEA}/genes/convergence/results/evidence_matrix.csv", low_memory=False)
    k = E.call_rate >= floor
    print(f"\ncandidate pool under a call-rate >= {floor} floor\n")
    t = pd.DataFrame({
        "variants": E.groupby("n_lines").size(),
        "kept": E[k].groupby("n_lines").size(),
    }).fillna(0).astype(int)
    t["pct_kept"] = (100 * t.kept / t.variants).round(0)
    print(t.to_string())
    print(f"\ntotal {len(E):,} -> {int(k.sum()):,} ({100*k.mean():.1f}%)")
    drop5 = E[(E.n_lines >= 4) & ~k]
    if len(drop5):
        print(f"\ntop candidates that would be dropped ({len(drop5)}):")
        print(drop5[["chrom", "pos", "target_gene", "gea_best_axis", "call_rate", "n_lines"]]
              .sort_values("n_lines", ascending=False).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["posthoc", "pool"], default="posthoc")
    ap.add_argument("--classes", default="nonsnp,sv,smallindel,snp")
    ap.add_argument("--floor", type=float, default=FLOOR)
    a = ap.parse_args()
    if a.mode == "pool":
        pool(a.floor)
    else:
        posthoc(a.classes.split(","), a.floor)
