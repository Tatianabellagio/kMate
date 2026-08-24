#!/usr/bin/env python
"""Standalone dot-grid of the significant (SV / small-indel / non-SNP) blocks across the
20 climate axes — LFMM, raw (before-WZA), clq0.9 TILING partition.

Lifted OUT of `raw_manhattan_snp_vs_{cls}_lfmm_tile.ipynb` so the dot plot can be
regenerated on its own: that notebook is slow because it renders 20 dense per-axis
mirror Manhattans first. Here we skip all Manhattan rendering and go straight to the grid.

Shows **ALL** clq0.9 tiling blocks whose LEAD (min-p) class record clears Bonferroni on
>=1 axis, and MARKS each (gene, axis) dot by whether the co-located SNP is also significant:
  * black-edged / solid  = non-SNP-only (no significant SNP in that block),
  * white-edged / faded   = SNP also significant in the same block.
(The earlier version showed only the non-SNP-only subset; now the SNP-also blocks are
visible too, per request.)

Additions kept from before:
  * denser rows (less vertical whitespace),
  * a right-hand panel of the variant's indel/SV **size** (|alt_len-ref_len|, bp),
    with rows ordered by that size (size per gene = its strongest-signal lead record).

Emits twin notebooks (same code, class swapped):
  newpeak_dotgrid_lfmm_nonsnp.ipynb  · _sv.ipynb  · _smallindel.ipynb
"""
import sys
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = "lfmm"


def build(cls, title_word):
    C = []
    md = lambda s: C.append(new_markdown_cell(s))
    co = lambda s: C.append(new_code_cell(s))

    md(f"""# Significant {title_word} blocks across the 20 climate axes — LFMM, raw, clq0.9 tiling

Standalone dot grid, lifted out of `raw_manhattan_snp_vs_{cls}_lfmm_tile.ipynb` (which is
slow — it renders 20 dense mirror Manhattans first). Shows **ALL** clq0.9 **tiling** blocks
whose LEAD (min-p) {cls} record clears per-class Bonferroni (0.05/n) on >=1 climate axis
(bio1-19 + pc1) — gen9, MAF>0.05, `wza_in_clq09_tile` LFMM p.

Each (gene, axis) dot is **marked by whether the co-located SNP is also significant**:
* **black-edged, solid** = **non-SNP-only** — no significant SNP in that block ("what did
  {title_word} find that SNPs missed"),
* **white-edged, faded** = **SNP also significant** in the same block (shown for context).

Dot **size** = -log10 p ({cls} lead); **color** = functional category; rows ordered by indel/SV
size (right panel). **Raw, uncalibrated p — inflation acknowledged, kept deliberately** (this is
the raw-LFMM working signal; see `nonsnp_block_characterization.ipynb` for the caveats).""")

    co(f"""import os, sys
import numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), "..", "..")))
import lib
GEA = lib.GEA
WZAIN = f"{{GEA}}/phase1_replication/results/multiaxis/wza_in_clq09_tile"
AXES = [f"bio{{i}}" for i in range(1, 20)] + ["pc1"]
MODEL = "{MODEL}"
CLS = "{cls}"
plt.rcParams.update({{"figure.dpi":110}})

def block_spans_tiling(r2=0.9):
    \"\"\"Tiling block spans: each block absorbs the preceding gap (matches reblock_blockdef.py).\"\"\"
    tag=f"clq{{r2}}"; rows=[]
    for ci in range(1,6):
        ch = f"Chr{{ci}}"
        g=pd.read_csv(f"{{lib.CLQ_BLOCKS_DIR}}/chr{{ci}}_{{tag}}_blocks_{{tag}}.tsv",sep="\\t").sort_values("start_pos").reset_index(drop=True)
        ends = g.end_pos.to_numpy(np.int64)
        for i in range(len(g)):
            lo = 1 if i == 0 else int(ends[i-1]) + 1
            hi = int(ends[i]) if i < len(g)-1 else 10**9
            rows.append((f"{{ch}}_{{i}}", ch, lo, hi))
    return pd.DataFrame(rows,columns=["block","chrom","start","end"]).set_index("block")
SPANS=block_spans_tiling(); GENES=lib.load_genes()
def genes_on(block):
    if block not in SPANS.index: return []
    s=SPANS.loc[block]; gc=GENES[(GENES.chrom==s.chrom)&(GENES.end>=s.start)&(GENES.start<=s.end)]
    return list(gc.gene)
print("loaded:", len(SPANS), "tiling blocks,", len(GENES), "genes ; class =", CLS)""")

    co("""# --- recompute all significant class blocks across all axes (no Manhattan rendering) ---
def load_raw(cls, axis):
    f = f"{WZAIN}/{MODEL}_{cls}_gen9_{axis}.csv"
    if not os.path.exists(f): return None
    d = pd.read_csv(f)
    d = d[d["MAF"] > 0.05].copy()
    d["nlp"] = -np.log10(d["pval"].clip(lower=1e-300))
    d["size"] = (d["alt_len"] - d["ref_len"]).abs()          # |indel/SV| length change, bp
    return d

def leads(d):
    d = d[d["block"].notna() & (d["block"] != "")]
    return d.loc[d.groupby("block")["nlp"].idxmax()]         # lead (min-p) record per block

def sig_blocks(axis, cls):
    s, c = load_raw("snp", axis), load_raw(cls, axis)
    if s is None or c is None:
        return None
    bonf_s, bonf_c = -np.log10(0.05 / len(s)), -np.log10(0.05 / len(c))
    ls_, lc_ = leads(s), leads(c)
    m = lc_.merge(ls_[["block", "nlp"]], on="block", suffixes=("", "_snp"), how="left")
    m["nlp_snp"] = m["nlp_snp"].fillna(0.0)                  # SNP-absent block -> lead SNP nlp 0
    sig = m[m["nlp"] > bonf_c].copy()                        # ALL class-significant blocks
    sig["snp_sig"] = sig["nlp_snp"] > bonf_s                 # is the co-located SNP also significant?
    sig["genes"] = sig["block"].map(lambda b: ";".join(genes_on(b)))
    sig["axis"] = axis
    return sig[["block", "chrom", "pos", "size", "nlp", "nlp_snp", "snp_sig", "genes", "axis"]]

all_sig = []
for a in AXES:
    x = sig_blocks(a, CLS)
    if x is None:
        print("missing", a); continue
    all_sig.append(x)
    n_only = int((~x["snp_sig"]).sum())
    print(f"{a}: {len(x)} sig {CLS} blocks ({n_only} non-SNP-only, {len(x)-n_only} SNP-also); "
          f"{(x.genes != '').sum()} gene-bearing")
ALLPK = pd.concat(all_sig, ignore_index=True) if all_sig else pd.DataFrame()
print(f"total (block,axis) significant: {len(ALLPK)} ; non-SNP-only {int((~ALLPK.snp_sig).sum())}")""")

    co("""# --- gene-level table: one row per (gene, axis), plus a representative (strongest) row per gene ---
GENE = ALLPK.assign(gene=ALLPK["genes"].str.split(";")).explode("gene")
GENE = GENE[GENE["gene"].astype(bool)].copy()
GA = GENE.sort_values("nlp", ascending=False).drop_duplicates(["gene", "axis"])   # 1 row / (gene,axis)
nax = GA.groupby("gene")["axis"].nunique().rename("n_axes")
nax_only = GA[~GA["snp_sig"]].groupby("gene")["axis"].nunique().rename("n_axes_nonsnp_only")
rep = (GA.sort_values("nlp", ascending=False).drop_duplicates("gene")             # strongest axis / gene
         .set_index("gene").join(nax).join(nax_only))
rep["n_axes_nonsnp_only"] = rep["n_axes_nonsnp_only"].fillna(0).astype(int)
print(f"{rep.shape[0]} unique gene-bearing significant {CLS} genes "
      f"({int((rep.n_axes_nonsnp_only>0).sum())} non-SNP-only on >=1 axis)")""")

    co('''# --- annotate genes (TAIR GO + UniProt) — the mandated annotator ---
import importlib.util as _ilu
_annp = os.path.abspath(os.path.join(os.getcwd(), "..", "annotate_genes_tair_uniprot.py"))
_sp = _ilu.spec_from_file_location("ann_tu", _annp)
ann_tu = _ilu.module_from_spec(_sp); _sp.loader.exec_module(ann_tu)

A = ann_tu.annotate(sorted(rep.index.tolist()))
A = A.set_index("gene")
rep = rep.join(A[["symbol", "protein_name", "categories"]])
rep["categories"] = rep["categories"].fillna("")
rep["symbol"] = rep["symbol"].fillna("")
_outp = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/newpeak_dotgrid_genes_{CLS}_tile.csv"
rep.reset_index().to_csv(_outp, index=False)
print(f"wrote {_outp}")
print("category counts (a gene can carry several):")
print(rep["categories"].str.split(",").explode().replace("", pd.NA).dropna().value_counts().to_string())''')

    co('''# --- dense dot grid, rows ordered by indel/SV size, with a right-hand size panel ---
import matplotlib.lines as mlines

AXES_ORD = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
CATCOL = {"flowering": "#8E44AD", "temperature": "#2980B9", "water": "#16A085",
          "light": "#F1C40F", "oxidative": "#E67E22", "defense": "#C0392B",
          "calcium": "#7F8C8D", "unclassified": "#CED4DA"}
def pcat(g):
    c = [x for x in (rep.loc[g, "categories"] or "").split(",") if x]
    return c[0] if c else "unclassified"

# show ALL gene-bearing significant genes (no recurrence subsampling). The figure is
# very tall for nonsnp/smallindel (~1,880 rows); ordered by size, NOT genomically.
keep = rep.index
gord = rep.loc[keep].sort_values(["size", "n_axes"]).index.tolist()          # small at bottom
thr = 1
yidx = {g: i for i, g in enumerate(gord)}
xidx = {a: i for i, a in enumerate(AXES_ORD)}
D = GA[GA["gene"].isin(keep)].drop_duplicates(["gene", "axis"])

H = max(4.0, 0.16 * len(gord))                                                # denser: 0.16 in/row
_dpi = min(110, int(62000 / H)) if H > 0 else 110                            # stay under Agg 65536px cap
fig = plt.figure(figsize=(11, H), dpi=_dpi)
gs = fig.add_gridspec(1, 2, width_ratios=[5.5, 1.0], wspace=0.04)
ax = fig.add_subplot(gs[0]); axr = fig.add_subplot(gs[1], sharey=ax)

for _, rr in D.iterrows():
    only = not rr["snp_sig"]                                  # non-SNP-only = no significant SNP in block
    ax.scatter(xidx[rr["axis"]], yidx[rr["gene"]], s=10 + rr["nlp"] * 2.0,
               c=CATCOL.get(pcat(rr["gene"]), "#CED4DA"),
               edgecolors=("black" if only else "white"),
               linewidths=(0.8 if only else 0.3),
               alpha=(0.98 if only else 0.45), zorder=(4 if only else 3))
ax.set_xticks(range(len(AXES_ORD))); ax.set_xticklabels(AXES_ORD, rotation=90, fontsize=7)
ax.set_yticks(range(len(gord)))
ax.set_yticklabels([f"{(rep.loc[g,'symbol'] or g)}  ({pcat(g)})" for g in gord], fontsize=6)
ax.set_xlim(-0.6, len(AXES_ORD) - 0.4); ax.set_ylim(-0.7, len(gord) - 0.3)
ax.grid(True, axis="both", lw=0.3, c="0.9", zorder=0); ax.set_axisbelow(True)
ax.set_xlabel("climate axis", fontsize=9)
for sp in ["top", "right"]:
    ax.spines[sp].set_visible(False)

# right panel: indel/SV size (bp) per gene, log x
sizes = [max(int(rep.loc[g, "size"]), 1) for g in gord]                       # floor 1 for log
axr.barh(range(len(gord)), sizes, height=0.72,
         color=[CATCOL.get(pcat(g), "#CED4DA") for g in gord], alpha=0.9, zorder=3)
axr.set_xscale("log")
axr.set_xlabel("|indel/SV| size (bp)", fontsize=9)
axr.grid(True, axis="x", lw=0.3, c="0.9", zorder=0); axr.set_axisbelow(True)
plt.setp(axr.get_yticklabels(), visible=False)
for sp in ["top", "right", "left"]:
    axr.spines[sp].set_visible(False)
axr.tick_params(axis="y", length=0)

_present = [c for c in CATCOL if c in {pcat(g) for g in gord}]
_h = [mlines.Line2D([], [], marker="o", ls="", ms=7, mfc=CATCOL[c], mec="white", label=c)
      for c in _present]
leg1 = ax.legend(handles=_h, title="category", loc="upper left", bbox_to_anchor=(1.28, 1.0),
                 frameon=False, fontsize=7, title_fontsize=8)
ax.add_artist(leg1)
_hm = [mlines.Line2D([], [], marker="o", ls="", ms=8, mfc="#9AA0A6", mec="black", mew=0.9,
                     label="non-SNP-only (no sig SNP)"),
       mlines.Line2D([], [], marker="o", ls="", ms=8, mfc="#9AA0A6", mec="white", mew=0.9,
                     alpha=0.5, label="SNP also significant")]
ax.legend(handles=_hm, title="block", loc="lower left", bbox_to_anchor=(1.28, 0.0),
          frameon=False, fontsize=7, title_fontsize=8)
fig.tight_layout(); plt.show()
_nonly = int((rep["n_axes_nonsnp_only"] > 0).reindex(gord).sum())
print(f"dot grid: {len(gord)} genes = ALL gene-bearing significant {CLS} blocks (no subsampling; "
      f"black-edged dot = non-SNP-only); {_nonly} of them non-SNP-only on >=1 axis.")''')

    nb = new_notebook(); nb["cells"] = C
    out = os.path.join(HERE, f"newpeak_dotgrid_lfmm_{cls}.ipynb")
    ep = ExecutePreprocessor(timeout=3600, kernel_name="basic", startup_timeout=180)
    ep.preprocess(nb, {"metadata": {"path": HERE}})
    with open(out, "w") as f:
        nbf.write(nb, f)
    print(f"[built+executed] {out}")


if __name__ == "__main__":
    todo = sys.argv[1:] or ["sv", "nonsnp"]
    names = {"sv": "SV", "nonsnp": "non-SNP (pooled)", "smallindel": "small-indel"}
    for cls in todo:
        build(cls, names.get(cls, cls))
