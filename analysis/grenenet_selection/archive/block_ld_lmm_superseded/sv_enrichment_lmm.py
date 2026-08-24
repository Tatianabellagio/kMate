#!/usr/bin/env python
"""SV enrichment in the CALIBRATED LD-LMM temporal-selection blocks (drift-null aware).

The magnitude-outlier version (site_sv_enrichment.py) called "selected" = top-1% by
|median SNP s|, which can include high-drift blocks. This uses instead the precision-
weighted LD-LMM (block_ld_lmm_temporal.py on the clq0.9 hapfreq matrix): per-block
selection with the linked background separated and a REML residual -> FDR q, lambda_GC
calibrated. "Selected" = LMM q<qcut. Same size-matched null (match n_nonsv) on the
same clq0.9 block grid -> is the answer (no SV enrichment) robust to the selection call?

INPUT: hapfreq_clq90/site{S}_block_ld_lmm_temporal.csv (LMM q per unit)
       site_temporal/site{S}_clq90_blocks.csv.gz (per-block class composition)
OUTPUT: site_temporal/site{S}_sv_enrichment_lmm.json
  PY=<plotting/basic>; $PY sv_enrichment_lmm.py --site 4 --qcut 0.05
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", type=int, default=4)
    ap.add_argument("--qcut", type=float, default=0.05)
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--min-snp", type=int, default=3)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    ST = f"{lib.GEA}/site_temporal"
    tab = pd.read_csv(f"{ST}/site{args.site}_clq90_blocks.csv.gz")
    lmm = pd.read_csv(f"{lib.GEA}/hapfreq_clq90/site{args.site}_block_ld_lmm_temporal.csv")

    # LMM-selected regions = q<qcut unit intervals; flag clq0.9 blocks overlapping them.
    sel_units = lmm[lmm.q < args.qcut][["chrom", "unit_start", "unit_end"]].drop_duplicates()
    selmask = np.zeros(len(tab), bool)
    posc = tab.set_index(np.arange(len(tab)))
    for c, g in sel_units.groupby("chrom"):
        bi = np.where(tab.chrom.values == c)[0]
        bs = tab.start_pos.values[bi]; be = tab.end_pos.values[bi]
        for us, ue in zip(g.unit_start.values, g.unit_end.values):
            selmask[bi[(bs <= ue) & (be >= us)]] = True
    tab["lmm_selected"] = selmask

    scored = tab[tab.n_snp >= args.min_snp].copy()
    sel = scored[scored.lmm_selected]; non = scored[~scored.lmm_selected]
    f_sel = float(sel.has_sv.mean()); f_non = float(non.has_sv.mean())
    n_sel = len(sel)

    rng = np.random.default_rng(args.seed)
    non_by = {nv: sub.index.to_numpy() for nv, sub in non.groupby("n_nonsv")}
    sizes = np.array(sorted(non_by))
    sel_sizes = sel.n_nonsv.to_numpy()
    null = np.empty(args.n_perm)
    for b in range(args.n_perm):
        picks = [rng.choice(non_by[sizes[np.argmin(np.abs(sizes - nv))]]) for nv in sel_sizes]
        null[b] = non.loc[picks].has_sv.mean()
    p_emp = float((np.sum(null >= f_sel) + 1) / (args.n_perm + 1))
    fold = f_sel / np.median(null)
    bal = float(non.loc[[rng.choice(non_by[sizes[np.argmin(np.abs(sizes - nv))]])
                         for nv in sel_sizes]].n_nonsv.mean())

    from scipy.stats import fisher_exact
    a, b_ = int(sel.has_sv.sum()), int((~sel.has_sv).sum())
    c, d = int(non.has_sv.sum()), int((~non.has_sv).sum())
    orr, praw = fisher_exact([[a, b_], [c, d]])

    print(f"site {args.site}: LMM q<{args.qcut} selected units={len(sel_units)} -> "
          f"clq0.9 blocks selected={n_sel} (scored {len(scored)})")
    print(f"  size balance n_nonsv: sel={sel.n_nonsv.mean():.1f} null={bal:.1f}")
    print(f"  has_SV  selected={f_sel:.3f}  null(med)={np.median(null):.3f}  "
          f"non-sel={f_non:.3f}  fold={fold:.2f}  p_emp={p_emp:.4f}")
    print(f"  raw 2x2 OR={orr:.2f} p={praw:.3g} ({a}/{a+b_} sel SV-bearing)")
    summ = dict(site=args.site, selection="calibrated LD-LMM clq0.9", qcut=args.qcut,
                lambda_gc=None, n_selected_units=int(len(sel_units)), n_selected_blocks=n_sel,
                hasSV_selected=f_sel, hasSV_null_median=float(np.median(null)),
                hasSV_nonselected=f_non, fold=float(fold), p_emp=p_emp,
                raw_OR=float(orr), raw_p=float(praw))
    json.dump(summ, open(f"{ST}/site{args.site}_sv_enrichment_lmm.json", "w"), indent=2)
    print(f"-> {ST}/site{args.site}_sv_enrichment_lmm.json")


if __name__ == "__main__":
    main()
