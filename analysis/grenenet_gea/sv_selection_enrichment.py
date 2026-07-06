#!/usr/bin/env python
"""Are SVs under selection? SV enrichment in block-specific-selected haploblocks (site 4).

For each dynld unit we count variants by class (SNP / SV / small-indel) by POSITION, mark units
as SELECTED if they contain an FDR (q<0.05) or Bonferroni (p<0.05/M) block-specific haplotype
from the temporal LD-LMM, and test whether SVs are over-represented in selected vs the tested
background -- CONTROLLING for block size (compare SV *fraction*, since bigger blocks carry more
of everything). Caveat: "block carries SV" != "SV is the causal allele" (within-block linkage);
this is the tractable block-level test of SV involvement in selection.

Env: kmate. SITE via env. Writes site<ID>_sv_enrichment.csv + figure + prints the verdict.
"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = "results/grenenet_gea/hapfreq"
CM = "results/grenenet_gea/phase1_replication/class_matrices"
SITE = int(os.environ.get("SITE", 4))
CLASSES = ["snp", "sv", "smallindel"]


def nchrom(c):
    return "Chr" + str(c).replace("Chr", "")


def main():
    reg = pd.read_csv(f"{H}/hapfreq_registry.csv")
    units = reg.drop_duplicates("unit_idx")[
        ["unit_idx", "chrom", "unit_start", "unit_end", "unit_nvar"]].copy().reset_index(drop=True)
    units["chrom"] = units.chrom.map(nchrom)
    N = len(units)

    # --- count variants of each class per unit, by position ---
    cnt = {c: np.zeros(N, int) for c in CLASSES}
    for c in CLASSES:
        rec = pd.read_csv(f"{CM}/{c}_gen9.records.csv", usecols=["chrom", "pos"])
        rec["chrom"] = rec.chrom.map(nchrom)
        for ch, gu in units.groupby("chrom"):
            ss = gu.unit_start.values; ee = gu.unit_end.values; ui = gu.index.values
            o = np.argsort(ss); ss, ee, ui = ss[o], ee[o], ui[o]
            pos = rec.loc[rec.chrom == ch, "pos"].values
            if not len(pos):
                continue
            j = np.searchsorted(ss, pos, side="right") - 1
            ok = (j >= 0) & (pos <= ee[np.clip(j, 0, len(ee) - 1)])
            np.add.at(cnt[c], ui[j[ok]], 1)
    units["n_snp"] = cnt["snp"]; units["n_sv"] = cnt["sv"]; units["n_indel"] = cnt["smallindel"]
    units["n_var"] = units.n_snp + units.n_sv + units.n_indel
    units["has_sv"] = units.n_sv > 0
    units["sv_frac"] = units.n_sv / np.clip(units.n_var, 1, None)
    units["key"] = units.chrom + ":" + units.unit_start.astype(str) + "-" + units.unit_end.astype(str)

    # --- selected units from the temporal LD-LMM (tested background = units with a tested hap) ---
    res = pd.read_csv(f"{H}/site{SITE}_block_ld_lmm_temporal.csv")
    res["chrom"] = res.chrom.map(nchrom)
    res["key"] = res.chrom + ":" + res.unit_start.astype(str) + "-" + res.unit_end.astype(str)
    M = len(res); bonf = 0.05 / M
    tested = set(res.key)
    fdr = set(res.loc[res.q < 0.05, "key"]); bonfk = set(res.loc[res.p < bonf, "key"])
    U = units[units.key.isin(tested)].copy()
    U["sel_fdr"] = U.key.isin(fdr); U["sel_bonf"] = U.key.isin(bonfk)
    U.to_csv(f"{H}/site{SITE}_sv_enrichment.csv", index=False)

    print(f"site {SITE}: {len(U):,} tested dynld units | SVs genome: {units.n_sv.sum():,}")
    print(f"  SV-carrying units overall: {100*U.has_sv.mean():.0f}%")

    def report(selcol, lab):
        sel = U[U[selcol]]; bg = U[~U[selcol]]
        print(f"\n=== {lab}: {len(sel)} selected vs {len(bg):,} background ===")
        # size confound check
        print(f"  block size (n_var) median: selected {sel.n_var.median():.0f}  vs  bg {bg.n_var.median():.0f}")
        # presence (raw)
        a = int(sel.has_sv.sum()); b = len(sel) - a; c = int(bg.has_sv.sum()); dd = len(bg) - c
        OR, p = stats.fisher_exact([[a, b], [c, dd]])
        print(f"  carries >=1 SV: selected {100*a/len(sel):.0f}% ({a}/{len(sel)})  vs  bg {100*c/len(bg):.0f}%"
              f"  | Fisher OR={OR:.2f}, p={p:.2g}")
        # size-controlled: SV fraction (normalized for block size)
        mw = stats.mannwhitneyu(sel.sv_frac, bg.sv_frac, alternative="two-sided")
        print(f"  SV fraction (size-normalized) median: selected {sel.sv_frac.median():.3f}  vs  "
              f"bg {bg.sv_frac.median():.3f}  | Mann-Whitney p={mw.pvalue:.2g}")
        # size-matched presence: restrict bg to similar-size blocks
        lo, hi = sel.n_var.quantile([.1, .9])
        bgm = bg[bg.n_var.between(lo, hi)]
        if len(bgm):
            cm = bgm.has_sv.mean()
            ORm, pm = stats.fisher_exact([[a, b], [int(bgm.has_sv.sum()), len(bgm) - int(bgm.has_sv.sum())]])
            print(f"  size-matched bg ({len(bgm)} units, n_var in [{lo:.0f},{hi:.0f}]): "
                  f"carries-SV {100*cm:.0f}%  | OR={ORm:.2f}, p={pm:.2g}")
        return sel, OR, p

    sel_fdr, OR_f, p_f = report("sel_fdr", "FDR q<0.05")
    sel_bonf, OR_b, p_b = report("sel_bonf", "Bonferroni")

    # which selected blocks carry SVs
    print("\n  FDR-selected blocks carrying SVs (top by SV count):")
    sv_blocks = sel_fdr[sel_fdr.has_sv].sort_values("n_sv", ascending=False)
    print(sv_blocks.head(10)[["chrom", "unit_start", "unit_end", "n_var", "n_sv", "sv_frac"]].to_string(index=False))

    # --- figure: SV fraction selected vs background (size-controlled view) ---
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for a_, col, lab in [(ax[0], "sel_fdr", "FDR q<0.05"), (ax[1], "sel_bonf", "Bonferroni")]:
        data = [U.loc[~U[col], "sv_frac"], U.loc[U[col], "sv_frac"]]
        a_.boxplot(data, labels=["background", "selected"], showfliers=False, widths=.5)
        a_.scatter(np.full(U[col].sum(), 2) + np.random.RandomState(0).uniform(-.08, .08, U[col].sum()),
                   U.loc[U[col], "sv_frac"], s=10, c="#c0392b", alpha=.5, zorder=3)
        a_.set_ylabel("SV fraction of block variants"); a_.set_title(lab, fontsize=10, loc="left")
        a_.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Site {SITE}: are SVs enriched in block-specific-selected haploblocks?", fontsize=11)
    fig.tight_layout(); fig.savefig(f"{H}/site{SITE}_sv_enrichment.png", dpi=150)
    print(f"\n[done] {H}/site{SITE}_sv_enrichment.csv + figure")


if __name__ == "__main__":
    main()
