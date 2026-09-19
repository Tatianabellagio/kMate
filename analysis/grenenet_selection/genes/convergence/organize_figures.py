#!/usr/bin/env python
"""One tidy figure tree for the candidate review, built from the review tables.

Files are HARD LINKS into the render dirs, so this tree costs no extra disk and a
re-render followed by a re-run picks the new file up.

Why this exists: figures accumulated in four places -- results/plots/screen (round-1 and
round-2 screen renders), results/plots/screen_ownaxis (an earlier round-1 own-axis
re-render), results/plots/loci (older dissection locus figures) and
results/plots/screen_review_small (downscaled reading copies) -- because every plotting
script wrote to its own directory, and the dissection scripts write to their own
results/plots/ before being copied. On top of that the GWAS re-render added `_GWAS`
names beside the originals. Nothing tells you which file is the CURRENT allele.

This script copies, for every reviewed gene, the figures that match the representative
variant in the current tables into

  genes/figures/       ONE flat directory -- every figure for every gene, no subfolders.
                       <sym>__grid.png / <sym>__locus.png / <sym>__atac.png, so a gene's
                       figures sort together. Grouping (which scan, review round, top-20
                       rank) lives in INDEX.csv columns, not in the directory layout.
    INDEX.csv          gene, sym, set, rank, verdict, grade, axis, store_row, grid, locus, atac
    README.md          provenance + how to read the three figure types

File names are `<sym>__grid.png` (per-garden allele-frequency trajectories),
`<sym>__locus.png` (Manhattan + founder panel + LD) and `<sym>__atac.png` (functional
tracks). PDFs are linked next to the PNGs when they exist. The source directories are left
untouched: they stay the raw render output.

env: kmate (pandas only). Re-runnable: the tree is rebuilt from scratch each time.
"""
from __future__ import annotations
import os
import re
import shutil
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = f"{HERE}/results"
GENES = os.path.dirname(HERE)                  # analysis/grenenet_selection/genes
# Every render dir under genes/, whichever layer produced it. dissection/ is included
# because ~400 genes have a garden-trajectory figure there and nowhere else; without it
# the tree silently omits them.
SRC = [f"{RES}/plots/screen", f"{RES}/plots/screen_ownaxis", f"{RES}/plots/loci",
       f"{GENES}/dissection/results/plots"]
# One tree for all three layers, at the genes/ root rather than buried in
# convergence/results/, so every figure for every gene is in one place.
DST = f"{GENES}/figures"
# listing each source dir once: these live on Lustre, where a few thousand stat() calls
# cost more than the whole copy
LISTING = {d: (set(os.listdir(d)) if os.path.isdir(d) else set()) for d in SRC}


def link(src: str, dst: str) -> None:
    """Hard-link (same filesystem, no second copy of ~1 GB of figures); copy if refused."""
    if os.path.exists(dst):
        os.remove(dst)
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)

TOP20 = ["FUS3", "GPX6", "AT2G30000", "AT4G13200", "SSL7", "CRK14", "SCPL34", "GLR1.3",
         "LPP1", "SEC10a", "OHP1", "MQJ16.10", "HIGLE", "AT1G30820", "AT4G13495",
         "BRCC36A", "APY3", "AT1G51480", "PAE5", "CAM5"]


def find(stems, kind):
    """Newest existing file for any of `stems`, preferring the later-rendered one."""
    pats = {"grid": ["{}_garden_trajectories_ownaxis", "{}_garden_trajectories"],
            "locus": ["{}_combined_gwas", "{}_combined"]}[kind]
    hits = []
    for d in SRC:
        for st in stems:
            for p in pats:
                for ext in (".png", ".pdf"):
                    name = f"{p.format(st)}{ext}"
                    if name in LISTING[d]:
                        hits.append(f"{d}/{name}")
    if not hits:
        return []
    png = [h for h in hits if h.endswith(".png")]
    if not png:
        return []
    best = max(png, key=os.path.getmtime)
    out = [best]
    pdf = best[:-4] + ".pdf"
    if os.path.basename(pdf) in LISTING[os.path.dirname(pdf)]:
        out.append(pdf)
    return out


TAKEN = {}          # basename -> the gene that owns it, so two syms cannot collide


def copy_for(sym, gene, gwas):
    """Link grid + locus for one gene into the flat tree; returns the names written.

    Everything lands in ONE directory named `<sym>__<kind>.<ext>`, so all of a gene's
    figures sort together and there is a single place to look. Where the same symbol is
    carried by two different genes (round-1 and round-2 tables overlap), the second one
    takes `<sym>_<gene>__<kind>` rather than silently overwriting the first -- in the old
    per-folder tree they were kept apart by the folder and that hazard did not exist.
    """
    raw = str(sym)
    s = raw.replace("/", "_").replace(" ", "_")
    # some renders kept the symbol verbatim, spaces and all ("KAS III_garden_...png")
    stems = [f"{s}_GWAS", s, raw] if gwas else [s, raw, f"{s}_{gene}", f"{s}_GEA"]
    written = {}
    for kind in ("grid", "locus"):
        for f in find(stems, kind):
            ext = os.path.splitext(f)[1]
            name = f"{s}__{kind}{ext}"
            if TAKEN.get(name, gene) != gene:
                name = f"{s}_{gene}__{kind}{ext}"
            TAKEN[name] = gene
            link(f, f"{DST}/{name}")
            written.setdefault(kind, name)
    return written


def main():
    R = pd.read_csv(f"{RES}/screen_visual_review_round2.csv")
    V = pd.read_csv(f"{RES}/screen_visual_review.csv")
    V = V[~V.get("audit_2026_09_17", pd.Series("", index=V.index)).fillna("")
           .str.startswith("DROPPED")]
    if os.path.isdir(DST):
        shutil.rmtree(DST)
    os.makedirs(DST, exist_ok=True)

    rows = []
    done = R[R.status == "done"]
    for r in done.itertuples():
        gw = r.set == "GWAS"
        w = copy_for(r.sym, r.target_gene, gw)
        rows.append(dict(gene=r.target_gene, sym=r.sym, set=r.set, rank="", verdict=r.verdict,
                         grade=r.grade, axis=r.axis, store_row=r.store_row,
                         grid=w.get("grid", ""), locus=w.get("locus", ""), atac=""))
    for r in V.itertuples():
        w = copy_for(r.sym, r.target_gene, False)
        rows.append(dict(gene=r.target_gene, sym=r.sym, set="round1", rank="",
                         verdict=r.final_tier, grade=r.own_verdict, axis=r.gea_best_axis,
                         store_row=r.store_row,
                         grid=w.get("grid", ""), locus=w.get("locus", ""), atac=""))
    I = pd.DataFrame(rows)
    for i, sym in enumerate(TOP20, 1):
        m = I[(I.sym == sym) & (I.set != "round1")]
        if not len(m):
            m = I[I.sym == sym]
        if not len(m):
            print(f"  top20: no figures for {sym}")
            continue
        I.loc[m.index[0], "rank"] = str(i)          # rank is a column now, not a folder
    # Every render on disk, in any source dir, for a gene the review tables do NOT carry:
    # genes withdrawn on audit (AT5G44220), ones that failed the LD-confirm (BT4), the
    # dissection-era ROBUST calls that predate the screen (EMB1241, GSH1), and the CARK
    # figures. Without this pass they exist on disk but appear nowhere in the tree, which
    # is the whole failure this script was written to end.
    seen = set(I.sym.astype(str))
    urows, done = [], set()
    stems = set()
    for d in SRC:
        for f in LISTING[d]:
            if f.endswith(".png"):
                st = re.sub(r"_GWAS$|_GEA$", "",
                            re.sub(r"_garden_trajectories.*|_combined.*", "", f))
                if st and st not in seen:
                    stems.add(st)
    # Genes outside the review tables get the SAME <sym>__grid / <sym>__locus names as
    # everyone else, picked by the same newest-render rule (find). They used to keep their
    # raw render names, so a gene rendered twice showed up as two differently named grids.
    linked = set()
    for st in sorted(stems):
        w = {}
        for kind in ("grid", "locus"):
            for f in find([st, f"{st}_GWAS", f"{st}_GEA"], kind):
                name = f"{st}__{kind}{os.path.splitext(f)[1]}"
                link(f, f"{DST}/{name}"); TAKEN[name] = st
                w.setdefault(kind, name)
                done.add(os.path.basename(f))
        if w:
            linked.add(st)
            urows.append(dict(gene="", sym=st, set="unreviewed", rank="", verdict="",
                              grade="", axis="", store_row="", grid=w.get("grid", ""),
                              locus=w.get("locus", ""), atac=""))
    # anything left that is not a grid/locus render (e.g. cark_sv_garden_dynamics) keeps
    # its original name
    for d in SRC:
        for f in sorted(LISTING[d]):
            if not f.endswith((".png", ".pdf")) or f in done:
                continue
            st = re.sub(r"_GWAS$|_GEA$", "", re.sub(r"_garden_trajectories.*|_combined.*", "", f))
            # skip only stems that actually got __grid/__locus links; the CARK figures
            # match neither pattern and would otherwise vanish from the tree
            if st in seen or st in linked or not st:
                continue
            link(f"{d}/{f}", f"{DST}/{f}"); done.add(f)
            if f.endswith(".png"):
                urows.append(dict(gene="", sym=st, set="unreviewed", rank="", verdict="",
                                  grade="", axis="", store_row="", grid="", locus=f, atac=""))
    I = pd.concat([I, pd.DataFrame(urows)], ignore_index=True)

    # expression figures (plot_expression.py): carriers vs non-carriers + lineage control
    expr_src = f"{RES}/plots/expr"
    if os.path.isdir(expr_src):
        erows = []
        for f in sorted(os.listdir(expr_src)):
            ext = os.path.splitext(f)[1]
            if ext not in (".png", ".pdf") or not f.startswith("expr_"):
                continue
            sym = f[len("expr_"):-len(ext)]
            link(f"{expr_src}/{f}", f"{DST}/{sym}__expr{ext}")
            if ext == ".png":
                erows.append(dict(gene="", sym=sym, set="expr", rank="", verdict="", grade="",
                                  axis="", store_row="", grid="", locus="", atac="",
                                  expr=f"{sym}__expr.png"))
        I = pd.concat([I, pd.DataFrame(erows)], ignore_index=True)

    # functional-track figures (plot_atac.py) -- not per-gene review renders, so they are
    # linked wholesale rather than matched to a representative variant.
    atac_src = f"{RES}/plots/atac"
    if os.path.isdir(atac_src):
        arows = []
        for f in sorted(os.listdir(atac_src)):
            ext = os.path.splitext(f)[1]
            if ext not in (".png", ".pdf"):
                continue
            sym = f[len("atac_locus_"):-len(ext)] if f.startswith("atac_locus_") else ""
            name = f"{sym}__atac{ext}" if sym else f
            link(f"{atac_src}/{f}", f"{DST}/{name}")
            if ext == ".png" and sym:
                arows.append(dict(gene="", sym=sym, set="atac", rank="", verdict="",
                                  grade="", axis="", store_row="", grid="", locus="",
                                  atac=name))
        # concat, never rebuild from `rows` -- the TOP20 loop above writes `rank` into I,
        # and those edits are not in `rows`.
        I = pd.concat([I, pd.DataFrame(arows)], ignore_index=True)

    I.sort_values(["set", "rank", "sym"]).to_csv(f"{DST}/INDEX.csv", index=False)

    files = os.listdir(DST)
    n = {"png": sum(f.endswith(".png") for f in files),
         "pdf": sum(f.endswith(".pdf") for f in files)}
    with open(f"{DST}/README.md", "w") as fh:
        fh.write(f"""# Candidate figures -- one flat folder

Every figure for every gene under `genes/`, in one place. Built by
`convergence/organize_figures.py`; re-run it after any re-render. Entries are hard links,
so this costs no extra disk and the raw render output stays where the scripts put it
(`convergence/results/plots/{{screen,screen_ownaxis,loci,atac}}`,
`dissection/results/plots`). Every reviewed file matches the representative variant in the
CURRENT tables (allele-resolved for the GWAS set, GEA-significant records only).

FLAT -- one directory, {n['png']} PNGs and {n['pdf']} PDFs, no subfolders. A gene's figures
share its symbol prefix, so they sort together and `ls | grep GPX6` finds all of them.
Grouping that used to be folders (top-20 rank, which scan, review round) is in the
`set` and `rank` columns of INDEX.csv instead.

Three figure types per gene:

  `<sym>__grid.png`   one panel per garden, ordered by the climate axis; dots = pools,
                      bold line = garden mean, dashed = founding frequency, red frame +
                      star = Bonferroni-significant in that garden's per-garden GWAS.
                      This is the allele-frequency-over-time view.
  `<sym>__locus.png`  top: -log10 p around the gene (lead variant outlined, gene shaded,
                      dashed lines = SNP / SV Bonferroni); middle: founder genotypes
                      sorted by home temperature; bottom: founder LD.
  `<sym>__atac.png`   functional tracks from `plot_atac.py`: GEA Manhattan, TFBS turnover,
                      per-tissue open chromatin, gene models on a shared axis.

INDEX.csv -- gene, sym, set, rank, verdict, grade, axis, store_row, grid, locus, atac.

    set=GEA_ownaxis  211 GEA genes judged on their own climate axis
    set=GWAS          45 per-garden GWAS genes
    set=round1        82 round-1 genes, the 7 withdrawn on 2026-09-17 removed
    set=unreviewed    rendered but carried by no review table: genes withdrawn on audit
                      (AT5G44220), ones that failed the LD-confirm (BT4), the
                      dissection-era ROBUST calls (EMB1241, GSH1), the CARK figures
    set=atac          functional-track figures
    rank=1..20        the genes to look at first

Files whose names are not `<sym>__*` are `set=unreviewed`: they keep their original render
names because they are not grid/locus review pairs. The coverage check this tree must pass
is that every gene with a PNG in any source dir has a row in INDEX.csv -- `unreviewed` is
what makes that true, without it ~90 genes existed on disk and nowhere here.

Verdicts: A pattern + co-signal locus, B pattern + lone SNP-blind variant, C pattern but
the lead is a passenger, D partial, E none. All are single-pass visual calls on raw
uncalibrated LFMM / GEMMA p-values.
""")
    print(f"wrote {DST} (flat): {n['png']} png + {n['pdf']} pdf")
    print(f"INDEX.csv rows: {len(I)}; genes with no grid found: "
          f"{int(((I.grid == '') & (~I.set.isin(['atac', 'unreviewed', 'expr']))).sum())}")


if __name__ == "__main__":
    main()
