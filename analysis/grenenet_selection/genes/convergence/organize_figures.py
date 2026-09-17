#!/usr/bin/env python
"""One tidy figure tree for the candidate review, built from the review tables.

Files are HARD LINKS into results/plots/, so this tree costs no extra disk and a
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

  results/figures/
    01_top20/            rank-ordered, the shortlist handed to the user
    gea_own_axis/        all 211 GEA own-axis genes (round 2)
    gwas/                all 45 GWAS genes (allele-resolved, `_GWAS` renders)
    round1/              the 89 round-1 genes still standing
    INDEX.csv            gene, set, rank, verdict, grade, axis, files, store_row
    README.md            provenance + how to read the two figure types

File names are `<sym>__grid.png` (per-garden trajectories) and `<sym>__locus.png`
(Manhattan + founder panel + LD). PDFs are copied next to the PNGs when they exist.
The source directories are left untouched: they stay the raw render output.

env: kmate (pandas only). Re-runnable: the tree is rebuilt from scratch each time.
"""
from __future__ import annotations
import os
import shutil
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = f"{HERE}/results"
SRC = [f"{RES}/plots/screen", f"{RES}/plots/screen_ownaxis", f"{RES}/plots/loci"]
DST = f"{RES}/figures"
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


def copy_for(sym, gene, gwas, folder, prefix=""):
    """Copy grid + locus for one gene; returns the names written."""
    raw = str(sym)
    s = raw.replace("/", "_").replace(" ", "_")
    # some renders kept the symbol verbatim, spaces and all ("KAS III_garden_...png")
    stems = [f"{s}_GWAS", s, raw] if gwas else [s, raw, f"{s}_{gene}", f"{s}_GEA"]
    written = {}
    for kind in ("grid", "locus"):
        for f in find(stems, kind):
            name = f"{prefix}{s}__{kind}{os.path.splitext(f)[1]}"
            link(f, f"{folder}/{name}")
            written.setdefault(kind, name)
    return written


def main():
    R = pd.read_csv(f"{RES}/screen_visual_review_round2.csv")
    V = pd.read_csv(f"{RES}/screen_visual_review.csv")
    V = V[~V.get("audit_2026_09_17", pd.Series("", index=V.index)).fillna("")
           .str.startswith("DROPPED")]
    if os.path.isdir(DST):
        shutil.rmtree(DST)
    for d in ("01_top20", "gea_own_axis", "gwas", "round1"):
        os.makedirs(f"{DST}/{d}", exist_ok=True)

    rows = []
    done = R[R.status == "done"]
    for r in done.itertuples():
        gw = r.set == "GWAS"
        sub = "gwas" if gw else "gea_own_axis"
        w = copy_for(r.sym, r.target_gene, gw, f"{DST}/{sub}")
        rows.append(dict(gene=r.target_gene, sym=r.sym, set=r.set, rank="", verdict=r.verdict,
                         grade=r.grade, axis=r.axis, store_row=r.store_row,
                         folder=sub, grid=w.get("grid", ""), locus=w.get("locus", "")))
    for r in V.itertuples():
        w = copy_for(r.sym, r.target_gene, False, f"{DST}/round1")
        rows.append(dict(gene=r.target_gene, sym=r.sym, set="round1", rank="",
                         verdict=r.final_tier, grade=r.own_verdict, axis=r.gea_best_axis,
                         store_row=r.store_row, folder="round1",
                         grid=w.get("grid", ""), locus=w.get("locus", "")))
    I = pd.DataFrame(rows)
    for i, sym in enumerate(TOP20, 1):
        m = I[(I.sym == sym) & (I.set != "round1")]
        if not len(m):
            m = I[I.sym == sym]
        if not len(m):
            print(f"  top20: no figures for {sym}")
            continue
        r = m.iloc[0]
        copy_for(sym, r.gene, r.set == "GWAS", f"{DST}/01_top20", prefix=f"{i:02d}_")
        I.loc[m.index[0], "rank"] = str(i)
    I.sort_values(["set", "rank", "sym"]).to_csv(f"{DST}/INDEX.csv", index=False)

    n = {d: len(os.listdir(f"{DST}/{d}")) for d in ("01_top20", "gea_own_axis", "gwas", "round1")}
    with open(f"{DST}/README.md", "w") as fh:
        fh.write(f"""# Candidate figures, one tree

Built by `organize_figures.py` from the review tables; re-run it after any re-render.
The raw render output stays in `../plots/` (screen/, screen_ownaxis/, loci/) -- this tree
is the curated view, and every file here matches the representative variant in the
CURRENT tables (allele-resolved for the GWAS set, GEA-significant records only).

    01_top20/       {n['01_top20']} files -- the 20 genes to look at first, rank-ordered
    gea_own_axis/   {n['gea_own_axis']} files -- 211 GEA genes judged on their own climate axis
    gwas/           {n['gwas']} files -- 45 per-garden GWAS genes
    round1/         {n['round1']} files -- the round-1 89, minus the 7 withdrawn on 2026-09-17
    INDEX.csv       gene, set, rank, verdict, grade, axis, store_row, file names

Two figure types per gene:

  `<sym>__grid.png`   one panel per garden, ordered by the climate axis; dots = pools,
                      bold line = garden mean, dashed = founding frequency, red frame +
                      star = Bonferroni-significant in that garden's per-garden GWAS.
  `<sym>__locus.png`  top: -log10 p around the gene (lead variant outlined, gene shaded,
                      dashed lines = SNP / SV Bonferroni); middle: founder genotypes
                      sorted by home temperature; bottom: founder LD.

Verdicts: A pattern + co-signal locus, B pattern + lone SNP-blind variant, C pattern but
the lead is a passenger, D partial, E none. All are single-pass visual calls on raw
uncalibrated LFMM / GEMMA p-values.
""")
    print(f"wrote {DST}: " + ", ".join(f"{k} {v} files" for k, v in n.items()))
    print(f"INDEX.csv rows: {len(I)}; genes with no grid found: "
          f"{int((I.grid == '').sum())}")


if __name__ == "__main__":
    main()
