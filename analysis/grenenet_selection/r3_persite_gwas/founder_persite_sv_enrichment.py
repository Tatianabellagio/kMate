#!/usr/bin/env python
"""Does the founder GWAS see SVs LOCALLY (per site) or only after cross-site integration?

Resolves the paradox: the founder-GWAS JOINT (SV x2.3) is a META across 30 sites. But it
KEEPS per-site EMMAX z (Z[block,site], kinship-corrected founder selection). So we can run the
IDENTICAL size-matched SV-enrichment three ways on the same clq0.9 blocks and compare:

  (a) founder-GWAS PER-SITE : selected = top-frac blocks by |Z[:,site]|  (per-site, kinship-corr)
  (b) founder-GWAS JOINT    : selected = top-frac blocks by chi2_joint   (cross-site meta)
  (c) [external] pool-temporal per-site  = cross_site_enrichment.csv (my |median SNP s|)

If (a) ~ (c) ~ 1.0 (n.s.) but (b) is enriched -> SVs emerge only from cross-site INTEGRATION,
not locally (answer: "selected at SOME site" borrowing strength, = local adaptation of DIFFERENT
SVs at different sites). If (a) >> (c) -> the kinship-corrected founder method sees per-site SV
selection that the crude per-variant |s| misses (answer: method, not integration).

INPUT: hapfreq/multisite_founder_gwas_clq90_pc1.{npz,csv}; site_temporal/site4_clq90_blocks.csv.gz
OUTPUT: site_temporal/founder_persite_sv_enrichment.csv (+ printed comparison)
  PY=<plotting>; $PY founder_persite_sv_enrichment.py --top-frac 0.02
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

ST = f"{lib.GEA}/site_temporal"
FG = f"{lib.GEA}/hapfreq/multisite_founder_gwas_clq90_pc1"


def matched_fold(scored, score, sel_mask, rng, n_perm):
    """Size-matched (n_nonsv) SV-bearing fold + empirical p for a selection mask."""
    sel = scored[sel_mask]; non = scored[~sel_mask]
    if len(sel) == 0 or len(non) == 0:
        return np.nan, np.nan, np.nan
    f_sel = sel.has_sv.mean()
    non_by = {nv: sub.index.to_numpy() for nv, sub in non.groupby("n_nonsv")}
    sizes = np.array(sorted(non_by)); hs = non.has_sv
    sel_sizes = sel.n_nonsv.to_numpy()
    null = np.empty(n_perm)
    for b in range(n_perm):
        picks = [rng.choice(non_by[sizes[np.argmin(np.abs(sizes - nv))]]) for nv in sel_sizes]
        null[b] = hs.loc[picks].mean()
    p = (np.sum(null >= f_sel) + 1) / (n_perm + 1)
    return float(f_sel), float(f_sel / np.median(null)), float(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top-frac", type=float, default=0.02)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--min-snp", type=int, default=3)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    z = np.load(f"{FG}.npz", allow_pickle=True)
    Z = z["Z"]; sites = z["sites"]; bio1 = z["bio1"]
    fg = pd.DataFrame(dict(chrom=z["chrom"].astype(str), start=z["start"], end=z["end"]))
    csv = pd.read_csv(f"{FG}.csv")
    # aggregate hap-clusters -> unit: per site max|z|, and max chi2_joint
    unit_key = fg.chrom + ":" + fg.start.astype(str) + "-" + fg.end.astype(str)
    absZ = np.abs(Z)
    df = pd.DataFrame(absZ, columns=[f"s{int(s)}" for s in sites])
    df["unit"] = unit_key.values; df["chi2_joint"] = csv["chi2_joint"].values
    agg = df.groupby("unit").max()  # max over clusters within a unit

    blk = pd.read_csv(f"{ST}/site4_clq90_blocks.csv.gz")
    blk["unit"] = blk.chrom + ":" + blk.start_pos.astype(str) + "-" + blk.end_pos.astype(str)
    B = agg.merge(blk[["unit", "n_nonsv", "n_snp", "has_sv"]], on="unit")
    scored = B[B.n_snp >= args.min_snp].reset_index(drop=True)
    print(f"scored clq0.9 units (>= {args.min_snp} SNP): {len(scored):,}  "
          f"has_sv rate {scored.has_sv.mean():.3f}")

    rng = np.random.default_rng(args.seed)
    # (b) JOINT
    k = max(1, int(round(args.top_frac * len(scored))))
    jmask = scored.chi2_joint.rank(ascending=False, method="first") <= k
    fj, foldj, pj = matched_fold(scored, scored.chi2_joint, jmask.to_numpy(), rng, args.n_perm)
    print(f"\n[JOINT meta]  top-{args.top_frac:.0%} ({k} units): "
          f"has_SV={fj:.3f} fold={foldj:.2f} p_emp={pj:.3f}")

    # (a) per-site
    rows = []
    for s in sites:
        col = f"s{int(s)}"
        mask = (scored[col].rank(ascending=False, method="first") <= k).to_numpy()
        f_sel, fold, p = matched_fold(scored, scored[col], mask, rng, args.n_perm)
        rows.append(dict(site=int(s), bio1=float(bio1[list(sites).index(s)]),
                         hasSV=f_sel, fold=fold, p_emp=p))
    R = pd.DataFrame(rows).sort_values("bio1")
    # merge my pool-temporal per-site fold
    try:
        pt = pd.read_csv(f"{ST}/cross_site_enrichment.csv")[["site", "fold", "p_emp"]]
        pt = pt.rename(columns={"fold": "pooltemporal_fold", "p_emp": "pooltemporal_p"})
        R = R.merge(pt, on="site", how="left")
    except Exception:
        pass
    R.to_csv(f"{ST}/founder_persite_sv_enrichment.csv", index=False)

    print(f"\n[founder-GWAS PER-SITE]  top-{args.top_frac:.0%}, size-matched SV-bearing fold:")
    print(R.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    n = len(R)
    print(f"\nfounder per-site: fold>1 {int((R.fold>1).sum())}/{n}, sig(p<.05) {int((R.p_emp<0.05).sum())}, "
          f"median fold {R.fold.median():.3f}")
    if "pooltemporal_fold" in R:
        pt = R.dropna(subset=["pooltemporal_fold"])
        print(f"pool-temporal per-site: fold>1 {int((pt.pooltemporal_fold>1).sum())}/{len(pt)}, "
              f"median {pt.pooltemporal_fold.median():.3f}")
        print(f"per-site founder vs pool-temporal fold corr (Spearman): "
              f"{pt[['fold','pooltemporal_fold']].corr(method='spearman').iloc[0,1]:.3f}")
    print(f"\n=> JOINT fold {foldj:.2f} (p={pj:.3f}) vs per-site median {R.fold.median():.3f} "
          f"=> {'INTEGRATION reveals SVs' if foldj>1.15 and R.fold.median()<1.1 else 'see interpretation'}")


if __name__ == "__main__":
    main()
