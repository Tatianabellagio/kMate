#!/usr/bin/env python
"""Consolidate the site-collapsed new-peak gene lists (SV + non-SNP, Bonferroni + FDR),
re-run them through the TAIR-GO/UniProt annotator, and summarise functional THEMES
(flowering/circadian, temperature, water, light, oxidative, defense, calcium) — cross-tabbed
against AXIS QUALITY, so an inflated-axis or pc1(structure) hit is never mistaken for signal.
"""
import os, sys, importlib.util as ilu
import numpy as np, pandas as pd
from scipy import stats

GEA = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MA = f"{GEA}/phase1_replication/results/multiaxis"
SITE = f"{MA}/wza_in_clq09_tile_site"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]

LISTS = {
    "sv_bonf":     f"{MA}/raw_manhattan_newpeak_genes_sv_site_tile_annotated.csv",
    "sv_fdr":      f"{MA}/raw_manhattan_newpeak_genes_sv_site_tile_FDR.csv",
    "nonsnp_bonf": f"{MA}/raw_manhattan_newpeak_genes_nonsnp_site_tile_annotated.csv",
    "nonsnp_fdr":  f"{MA}/raw_manhattan_newpeak_genes_nonsnp_site_tile_FDR.csv",
}
THEMES = ["flowering", "temperature", "water", "light", "oxidative", "defense", "calcium"]


def lam(p):
    p = np.clip(np.asarray(p, float), 1e-300, 1.0)
    return float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.isf(0.5, 1))


def axis_lambda(cls):
    out = {}
    for a in AXES:
        f = f"{SITE}/lfmm_{cls}_gen9_{a}.csv"
        if os.path.exists(f):
            d = pd.read_csv(f); out[a] = lam(d[d.MAF > 0.05].pval)
    return out


def quality(cls_axis_lam, axis):
    if axis == "pc1":
        return "pc1(structure)"
    L = cls_axis_lam.get(axis, np.nan)
    if not np.isfinite(L):
        return "?"
    return "clean" if L < 1.6 else ("watch" if L < 2.0 else "INFLATED")


def main():
    frames = []
    for name, f in LISTS.items():
        if not os.path.exists(f):
            print(f"[skip missing] {name}: {os.path.basename(f)}"); continue
        d = pd.read_csv(f); d["list"] = name
        d["cls"] = "sv" if name.startswith("sv") else "nonsnp"
        d["sig"] = d["q"] if "q" in d.columns else np.nan
        frames.append(d[["gene", "block", "chrom", "axis", "n_axes", "list", "cls", "sig"]
                        + [c for c in ["nlp"] if c in d.columns]])
    if not frames:
        print("no lists found yet"); return
    G = pd.concat(frames, ignore_index=True)
    lam_sv, lam_ns = axis_lambda("sv"), axis_lambda("nonsnp")
    G["axis_qual"] = [quality(lam_sv if c == "sv" else lam_ns, a) for c, a in zip(G.cls, G.axis)]

    union = sorted(G.gene.dropna().unique())
    print(f"loaded {len(G)} (gene,list,axis) rows | {len(union)} unique genes | annotating via TAIR/UniProt ...",
          flush=True)
    sp = ilu.spec_from_file_location("ann", f"{GEA}/phase1_replication/annotate_genes_tair_uniprot.py")
    ann = ilu.module_from_spec(sp); sp.loader.exec_module(ann)
    A = ann.annotate(union)
    keep = [c for c in ["gene", "protein_name", "categories", "climate_stress_flowering",
                        "uniprot_function", "uniprot_keywords"] if c in A.columns]
    A = A[keep].copy(); A["categories"] = A["categories"].fillna("")
    M = G.merge(A, on="gene", how="left")
    M.to_csv(f"{os.path.dirname(os.path.abspath(__file__))}/functional_themes_site.csv", index=False)

    # ---- theme x axis-quality cross-tab (unique genes) ----
    gene_cat = A.set_index("gene")["categories"].to_dict()
    gq = M.groupby("gene")["axis_qual"].agg(lambda s: "clean" if "clean" in set(s) else
                                            ("watch" if "watch" in set(s) else sorted(set(s))[0]))
    print("\n================  THEME x BEST-AXIS-QUALITY  (unique genes)  ================")
    hdr = f"{'theme':12s} {'clean':>6s} {'watch':>6s} {'INFLATED':>9s} {'pc1':>5s}"
    print(hdr)
    for th in THEMES:
        gs = [g for g in union if th in (gene_cat.get(g, "") or "")]
        cnt = {"clean": 0, "watch": 0, "INFLATED": 0, "pc1(structure)": 0}
        for g in gs:
            cnt[gq.get(g, "?")] = cnt.get(gq.get(g, "?"), 0) + 1
        print(f"{th:12s} {cnt['clean']:6d} {cnt['watch']:6d} {cnt['INFLATED']:9d} {cnt['pc1(structure)']:5d}")

    # ---- the genes worth taking seriously: a climate theme on a clean/watch axis ----
    print("\n================  THEMED genes on CLEAN or WATCH axes (not INFLATED, not pc1)  ================")
    good = M[M.axis_qual.isin(["clean", "watch"]) & (M.categories.str.len() > 0)].copy()
    good = good.sort_values(["axis_qual", "sig"], ascending=[True, True])
    cols = ["gene", "protein_name", "categories", "cls", "list", "axis", "axis_qual", "n_axes", "sig"]
    if len(good):
        print(good[cols].drop_duplicates(["gene", "axis"]).to_string(index=False))
    else:
        print("  (none — every themed gene sits on an INFLATED axis or pc1)")
    print(f"\nwrote functional_themes_site.csv ({len(union)} genes)")


if __name__ == "__main__":
    main()
