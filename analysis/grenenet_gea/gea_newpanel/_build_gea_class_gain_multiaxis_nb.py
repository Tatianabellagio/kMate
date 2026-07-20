#!/usr/bin/env python
"""Build+execute notebooks/gea_class_gain_multiaxis.ipynb (run in `basic` env), AFTER
the multiaxis_fresh array (bio2..bio19) completes. Aggregates the per-axis clq0.9-block
WZA climate-GEA into the 'genes gained from non-SNP / SV' tables, ONE ROW PER (axis, model),
for all bioclim axes bio1..bio19, all 3 models, 3 classes. Then the union gained-gene list
across all axes, annotated with Ensembl-Plants functions.

WZA inputs:
  bio1  : results/.../gea_newpanel/wza_{model}_{cls}_clq09.csv            (already fresh)
  bio2-19: results/.../gea_newpanel/multiaxis_fresh/wza_{model}_{cls}_{axis}_clq09.csv
"""
import os, sys, subprocess
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel"
OUT = f"{HERE}/notebooks/gea_class_gain_multiaxis.ipynb"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

md_title = r"""# Climate-GEA gained-peaks / gained-genes — ALL bioclim axes (bio1..bio19)

Per (climate axis, model) row: WZA-corrected clq0.9-block significant peaks by class
(SNP / non-SNP / SV) and the gain over SNPs (`n_nonsnp_not_snp`, `n_sv_not_snp` + pct),
at Bonferroni and BH-FDR. Then the union list of genes under non-SNP-only / SV-only peaks
across all axes, annotated with Ensembl-Plants function. LFMM is the calibrated model;
binomial/Kendall are inflated references (read as upper bounds)."""

code = r"""
import os, sys
import numpy as np, pandas as pd
os.chdir("/global/scratch/users/tbellg/kmate"); sys.path.insert(0, "analysis/grenenet_gea")
sys.path.insert(0, "analysis/grenenet_gea/gea_newpanel")
import lib, blocks_clq09
GNP = "analysis/grenenet_gea/gea_newpanel/results"; MA = f"{GNP}/multiaxis_fresh"
MODELS = ["binomial", "kendall", "lfmm"]; CLASSES = ["snp", "nonsnp", "sv"]
AXES = [f"bio{i}" for i in range(1, 20)]

def wza_path(model, cls, axis):
    return f"{GNP}/wza_{model}_{cls}_clq09.csv" if axis == "bio1" \
        else f"{MA}/wza_{model}_{cls}_{axis}_clq09.csv"

def sig_sets(model, cls, axis, q=0.05):
    f = wza_path(model, cls, axis)
    if not os.path.exists(f):
        return None
    w = pd.read_csv(f).rename(columns={"index": "block"}); w["block"] = w.block.astype(str)
    w = w[w.Z_pVal.notna()]
    if not len(w):
        return set(), set()
    p = w.Z_pVal.to_numpy()
    return set(w.loc[p < 0.05/len(w), "block"]), set(w.loc[lib.bh(p) < q, "block"])

# cache sig sets per (axis,model,class,threshold)
S = {}
avail = []
for ax in AXES:
    ok = all(sig_sets(m, c, ax) is not None for m in MODELS for c in CLASSES)
    if not ok:
        continue
    avail.append(ax)
    for m in MODELS:
        for c in CLASSES:
            b, f = sig_sets(m, c, ax); S[(ax, m, c, 0)] = b; S[(ax, m, c, 1)] = f
print(f"axes available: {len(avail)}/{len(AXES)} -> {avail}")
"""

code_tables = r"""
def gain_rows(thr):
    rows = []
    for ax in avail:
        for m in MODELS:
            s = {c: S[(ax, m, c, thr)] for c in CLASSES}
            row = {"axis": ax, "model": m, "n_snp": len(s["snp"]), "n_nonsnp": len(s["nonsnp"]), "n_sv": len(s["sv"])}
            for b in ("nonsnp", "sv"):
                new = s[b] - s["snp"]; row[f"n_{b}_not_snp"] = len(new)
                row[f"pct_{b}_not_snp"] = round(len(new)/len(s[b]), 3) if s[b] else np.nan
            rows.append(row)
    return pd.DataFrame(rows)

gain_bonf = gain_rows(0); gain_fdr = gain_rows(1)
gain_bonf.to_csv(f"{GNP}/snp_vs_nonsnp/gea_classgain_multiaxis_bonf.csv", index=False)
gain_fdr.to_csv(f"{GNP}/snp_vs_nonsnp/gea_classgain_multiaxis_fdr.csv", index=False)
print("per-(axis,model) tables:", gain_bonf.shape)
"""

code_genes = r"""
bmap = {r.block: (r.chrom, int(r.start_pos), int(r.end_pos)) for r in blocks_clq09.load_blocks().itertuples()}
genes = lib.load_genes(); FLANK = 2000

def gained_genes(thr):
    hits = {}
    for cls in ("nonsnp", "sv"):
        blk2 = {}
        for ax in avail:
            for m in MODELS:
                for b in (S[(ax, m, cls, thr)] - S[(ax, m, "snp", thr)]):
                    d = blk2.setdefault(b, {"axes": set(), "models": set()})
                    d["axes"].add(ax); d["models"].add(m)
        for b, info in blk2.items():
            if b not in bmap:
                continue
            c, s0, e0 = bmap[b]
            inb = genes[(genes.chrom == c) & (genes.start <= e0) & (genes.end >= s0)]
            flk = genes[(genes.chrom == c) & (genes.start <= e0+FLANK) & (genes.end >= s0-FLANK)]
            flk = flk[~flk.gene.isin(inb.gene)]
            for gid, otype in [(g, "in_block") for g in inb.gene] + [(g, "flank2kb") for g in flk.gene]:
                d = hits.setdefault((gid, cls), dict(gene=gid, klass=cls, chrom=c, overlap=otype,
                                                     blocks=set(), axes=set(), models=set()))
                d["blocks"].add(b); d["axes"].update(info["axes"]); d["models"].update(info["models"])
                if otype == "in_block":
                    d["overlap"] = "in_block"
    rows = [dict(gene=d["gene"], klass=d["klass"], chrom=d["chrom"], overlap=d["overlap"],
                 n_blocks=len(d["blocks"]), n_axes=len(d["axes"]), n_models=len(d["models"]),
                 axes=";".join(sorted(d["axes"])), models=";".join(sorted(d["models"]))) for d in hits.values()]
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["klass", "n_axes", "n_models", "gene"], ascending=[True, False, False, True]).reset_index(drop=True)
    return df

genes_bonf = gained_genes(0); genes_fdr = gained_genes(1)
gb = f"{GNP}/snp_vs_nonsnp/gea_gained_genes_multiaxis_bonf.csv"
gf = f"{GNP}/snp_vs_nonsnp/gea_gained_genes_multiaxis_fdr.csv"
genes_bonf.to_csv(gb, index=False); genes_fdr.to_csv(gf, index=False)
print(f"gained genes (union all axes) — Bonf {len(genes_bonf)}, FDR {len(genes_fdr)}")
# Ensembl-Plants function annotation (adds symbol + function columns in place)
import subprocess
subprocess.run([sys.executable, "analysis/grenenet_gea/gea_newpanel/annotate_gained_genes.py", gb, gf], check=False)
genes_bonf = pd.read_csv(gb).fillna(""); genes_fdr = pd.read_csv(gf).fillna("")
"""

nb = new_notebook(cells=[
    new_markdown_cell(md_title), new_code_cell(code),
    new_markdown_cell("## Table 1 — Bonferroni, per (axis, model)"), new_code_cell(code_tables + "\ngain_bonf"),
    new_markdown_cell("## Table 2 — BH-FDR, per (axis, model)"), new_code_cell("gain_fdr"),
    new_markdown_cell("## Gained genes across all axes (Ensembl-annotated)"), new_code_cell(code_genes),
    new_markdown_cell("### Genes — Bonferroni (non-SNP-only / SV-only, any axis)"), new_code_cell("genes_bonf"),
    new_markdown_cell("### Genes — FDR (non-SNP-only / SV-only, any axis)"), new_code_cell("genes_fdr"),
])
ep = ExecutePreprocessor(timeout=1800, kernel_name="python3", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": HERE}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
