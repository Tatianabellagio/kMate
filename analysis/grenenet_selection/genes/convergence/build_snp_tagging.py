#!/usr/bin/env python
"""Is each non-SNP candidate visible to SNPs? Two different questions, both answered.

The pool keeps SNP-shadowed and SNP-unique candidates alike (user, 2026-08-26: "only
hits that are non-SNPs, but including the ones that are shared with SNPs or not").
This script adds the columns that let you split them afterwards. They are NOT filters.

Two questions that are easy to conflate and must not be:

1. **Is the association signal also seen in SNPs?** -- `snp_cosig_2kb`: a genome-wide
   Bonferroni-significant SNP within 2 kb, on any of the 22 LFMM axes. This is
   statistical shadowing. If yes, the locus would have been found by a SNP study too;
   the non-SNP record is a better-resolved view of a known peak, not a new one.
   Method follows `dissection/augment_snp_cosig.py`, extended 20 -> 22 axes and with
   the Bonferroni threshold recomputed per axis instead of hardcoded.

2. **Could a SNP tag this variant at all?** -- `best_r2_snp`: the maximum founder-panel
   r2 between the variant and any SNP within +/-50 kb. This is genetic detectability,
   independent of whether anything was significant. A variant with low r2 to every
   nearby SNP is invisible to a SNP array by construction -- the "SNP-blind" layer that
   motivates the whole non-SNP track.

   Taken from the **masked** tagging run (`sv_snp_ld_v2/tagging_{sv,indel}_panel_*`),
   NOT `varexp/nonsnp_tagging_*`: the earlier build coded missingness as REF and had no
   MAC floor, and panel SV records average ~60% missing, so it inflated r2 for markers
   jointly missing in the same founders. `best_r2` is NaN where no SNP was testable in
   the window -- reported as `untestable`, not as "untagged", since the two mean
   opposite things.

A candidate can be any combination: SNP-shadowed and well-tagged (a known peak),
SNP-unique but well-tagged (a SNP array could see the locus but no SNP reached
significance), or SNP-unique and SNP-blind (the strongest case for the non-SNP layer
adding something).

The join is on **chrom + pos + size**, not chrom:pos. 2.14% of arch3 positions carry
several biallelic records (`panel-multiallelic-pos-key-trap`); adding size resolves
most of that ambiguity, and `tag_join_ambiguous` flags what it cannot.

Outputs: adds columns to results/gea_pool.csv and results/functional_variants.csv
(rewritten in place), and prints the shadowed/unique x tagged/blind breakdown.

env: kmate.  Run on a compute node.
"""
from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                    # noqa: E402
import axis_clusters as ac                                    # noqa: E402

WZAIN = (f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/"
         "wza_in_clq09_tile")
TAGDIR = f"{lib.GEA}/r3_persite_gwas/results/sv_snp_ld_v2"
OUT = f"{HERE}/results"

WIN = 2000              # SNP co-significance window, matches augment_snp_cosig.py
BLIND_R2 = 0.2          # "SNPs genuinely cannot see it", per build_nonsnp_tagging.py


def sig_snps() -> pd.DataFrame:
    """Every genome-wide Bonferroni SNP hit, across all 22 axes."""
    rows = []
    for ax in ac.AXES:
        f = f"{WZAIN}/lfmm_snp_gen9_{ax}.csv"
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f, usecols=["chrom", "pos", "MAF", "pval"])
        d = d[d.MAF > 0.05]
        thr = 0.05 / len(d)                    # per-axis, not a hardcoded constant
        s = d[d.pval < thr]
        if len(s):
            rows.append(pd.DataFrame({"chrom": s.chrom.to_numpy(str),
                                      "pos": s.pos.to_numpy(np.int64),
                                      "axis": ax}))
        print(f"  snp {ax:5s}  bonf<{thr:.3e}  sig={len(s):>5,}", flush=True)
    return pd.concat(rows, ignore_index=True)


def add_cosig(V: pd.DataFrame, S: pd.DataFrame) -> pd.DataFrame:
    """Nearest Bonferroni SNP within WIN bp, and which axes it is significant on."""
    flag, dist, axes = [], [], []
    by_chrom = {ch: d.sort_values("pos") for ch, d in S.groupby("chrom")}
    for r in V.itertuples():
        d = by_chrom.get(r.chrom)
        if d is None:
            flag.append(False); dist.append(np.nan); axes.append("")
            continue
        p = d.pos.to_numpy()
        lo, hi = int(r.pos) - WIN, int(r.pos) + WIN
        i, j = np.searchsorted(p, lo), np.searchsorted(p, hi, side="right")
        if j <= i:
            flag.append(False); dist.append(np.nan); axes.append("")
            continue
        sub = d.iloc[i:j]
        dd = np.abs(sub.pos.to_numpy() - int(r.pos))
        flag.append(True)
        dist.append(int(dd.min()))
        axes.append(",".join(sorted(set(sub.axis[dd == dd.min()]))))
    V = V.copy()
    V["snp_cosig_2kb"] = flag
    V["snp_cosig_dist"] = dist
    V["snp_cosig_axes"] = axes
    return V


def tagging_table() -> pd.DataFrame:
    """Masked founder-panel SV/indel -> SNP tagging r2, both classes, all chroms."""
    frames = []
    for kind in ("sv", "indel"):
        for ch in lib.CHROMS:
            f = f"{TAGDIR}/tagging_{kind}_panel_{ch}.npz"
            if not os.path.exists(f):
                print(f"  MISSING {os.path.basename(f)}", flush=True)
                continue
            d = np.load(f, allow_pickle=True)
            frames.append(pd.DataFrame({
                "chrom": d["chrom"].astype(str), "pos": d["pos"].astype(np.int64),
                "size": d["size"].astype(np.int64),
                "best_r2_snp": d["best_r2"], "best_r2_snp_pos": d["best_snp_pos"],
                "tag_n_snp_window": d["n_snp_window"],
                "tag_an": d["an"], "tag_maf": d["maf"]}))
    T = pd.concat(frames, ignore_index=True)
    # keep the best-tagged record per (chrom,pos,size); count how many collapsed
    T = T.sort_values("best_r2_snp", ascending=False)
    n = T.groupby(["chrom", "pos", "size"]).size().rename("tag_n_records")
    T = T.groupby(["chrom", "pos", "size"]).first().join(n).reset_index()
    T["tag_join_ambiguous"] = T.tag_n_records > 1
    return T


def annotate(path: str, S: pd.DataFrame, T: pd.DataFrame) -> pd.DataFrame:
    V = pd.read_csv(path)
    drop = [c for c in V.columns if c.startswith(("snp_cosig", "best_r2_snp", "tag_"))]
    V = V.drop(columns=drop)
    V = add_cosig(V, S)
    V = V.merge(T, on=["chrom", "pos", "size"], how="left")
    V["snp_blind"] = V.best_r2_snp < BLIND_R2
    V["tag_untestable"] = V.best_r2_snp.isna()
    V.to_csv(path, index=False)
    return V


def main():
    print("collecting Bonferroni SNP hits, 22 axes ...", flush=True)
    S = sig_snps()
    print(f"  {len(S):,} SNP hit-records, "
          f"{S.groupby(['chrom','pos']).ngroups:,} unique SNPs\n")

    print("loading masked panel tagging r2 ...", flush=True)
    T = tagging_table()
    print(f"  {len(T):,} tagged non-SNP records "
          f"({int(T.tag_join_ambiguous.sum()):,} ambiguous on chrom+pos+size)\n")

    for name in ("gea_pool.csv", "functional_variants.csv",
                 "functional_shortlist.csv"):
        p = f"{OUT}/{name}"
        if not os.path.exists(p):
            print(f"SKIP {name} (not built)")
            continue
        V = annotate(p, S, T)
        print(f"=== {name}: {len(V):,} rows ===")
        print(f"  SNP-shadowed (Bonferroni SNP within {WIN} bp): "
              f"{int(V.snp_cosig_2kb.sum()):,}  |  SNP-unique: "
              f"{int((~V.snp_cosig_2kb).sum()):,}")
        print(f"  tagging r2 recovered for {int(V.best_r2_snp.notna().sum()):,}; "
              f"untestable {int(V.tag_untestable.sum()):,}")
        tag = V[V.best_r2_snp.notna()]
        if len(tag):
            print(f"  median best_r2 {tag.best_r2_snp.median():.3f}; "
                  f"SNP-blind (r2<{BLIND_R2}): {int(tag.snp_blind.sum()):,}")
            print(pd.crosstab(tag.snp_cosig_2kb, tag.snp_blind,
                              rownames=["snp_cosig_2kb"],
                              colnames=["snp_blind"]).to_string())
        print()
    print(f"columns added in place under {OUT}/")


if __name__ == "__main__":
    main()
