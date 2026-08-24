#!/usr/bin/env python
"""Build+execute notebooks/gea_class_gain_tables.ipynb (run in `basic` env).

The 'genes gained from non-SNP / SV' tables, but on the CLIMATE GEA (not the per-garden
founder-selection GWAS): for the clq0.9-block WZA climate-GEA (gen9, bio1) it reports,
per model (binomial/kendall/lfmm), how many blocks each variant class (SNP / non-SNP /
SV) calls significant, and — the informative part — how many the non-SNP or SV scan
flags that the SNP scan does NOT (n_nonsnp_not_snp, n_sv_not_snp), at both Bonferroni
and BH-FDR. Then the genes under those non-SNP-only / SV-only blocks.

Same clq0.9 partition (blocks_recompute) + WZA (deg7-cap2000) as grid_3model; blocks
share ids across classes, so 'sig in class X not SNP' is a real spatial gain.

Inputs: analysis/grenenet_gea/gea_newpanel/results/wza_{model}_{cls}_clq09.csv  (model in
binomial/kendall/lfmm; cls in snp/nonsnp/sv). SV added by run_wza_sv_clq09.py.

  BPY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
  $BPY analysis/grenenet_gea/gea_newpanel/_build_gea_class_gain_tables_nb.py
"""
import os, sys
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel"
OUT = f"{HERE}/notebooks/gea_class_gain_tables.ipynb"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

md_title = r"""# Climate-GEA — SNP vs non-SNP vs SV: significant peaks & genes GAINED per layer

clq0.9-block WZA climate-GEA (gen9, **bio1**), three models (binomial / Kendall / LFMM)
x three variant classes (SNP / non-SNP = SV+indel / SV-only). Blocks share one clq0.9
partition across classes, so a block significant in non-SNP or SV but **not** in the SNP
scan is a genuine spatial gain from that layer.

- **Table 1 / Table 2** — per model: `n_snp`/`n_nonsnp`/`n_sv` significant blocks, and the
  gain columns `n_nonsnp_not_snp` / `n_sv_not_snp` (+ `pct` = fraction of that class's own
  hits that SNPs miss), at Bonferroni (0.05/n) and BH-FDR (q<0.05).
- **Gene tables** — the genes under non-SNP-only / SV-only significant blocks (union across
  models), annotated, with which models flagged them.

Reading guide (established in this project): binomial & Kendall are the **inflated raw
references** (pool pseudoreplication); **LFMM is the calibrated model**. Weight the LFMM
row / LFMM-flagged genes accordingly."""

code_load = r"""
import os, sys
import numpy as np, pandas as pd
os.chdir("/global/scratch/users/tbellg/kmate"); sys.path.insert(0, "analysis/grenenet_gea")
sys.path.insert(0, "analysis/grenenet_gea/gea_newpanel")
import lib, blocks_clq09
GNP = "analysis/grenenet_gea/gea_newpanel/results"
MODELS = ["binomial", "kendall", "lfmm"]; CLASSES = ["snp", "nonsnp", "sv"]

def load_wza(model, cls):
    w = pd.read_csv(f"{GNP}/wza_{model}_{cls}_clq09.csv").rename(columns={"index": "block"})
    w["block"] = w["block"].astype(str)
    return w[w["Z_pVal"].notna()].copy()

def sig_sets(model, cls, q=0.05):
    w = load_wza(model, cls); p = w["Z_pVal"].to_numpy()
    bonf = set(w.loc[p < 0.05 / len(w), "block"])
    qv = lib.bh(p)
    fdr = set(w.loc[qv < q, "block"])
    return bonf, fdr

# per (model, class) -> (bonf_set, fdr_set)
S = {(m, c): sig_sets(m, c) for m in MODELS for c in CLASSES}
print("WZA blocks per class:", {c: len(load_wza("lfmm", c)) for c in CLASSES})
"""

md_tables = r"""## Tables 1 & 2 — significant blocks per class + gain from non-SNP / SV"""

code_tables = r"""
def gain_table(thr):   # thr: 0=Bonferroni, 1=FDR
    rows = []
    for m in MODELS:
        s = {c: S[(m, c)][thr] for c in CLASSES}
        row = {"model": m, "n_snp": len(s["snp"]), "n_nonsnp": len(s["nonsnp"]), "n_sv": len(s["sv"])}
        for b in ("nonsnp", "sv"):
            new = s[b] - s["snp"]
            row[f"n_{b}_not_snp"] = len(new)
            row[f"pct_{b}_not_snp"] = round(len(new) / len(s[b]), 3) if s[b] else np.nan
        rows.append(row)
    return pd.DataFrame(rows)

gea_bonf_df = gain_table(0)
gea_fdr_df = gain_table(1)
gea_bonf_df.to_csv(f"{GNP}/snp_vs_nonsnp/gea_classgain_bonf.csv", index=False)
gea_fdr_df.to_csv(f"{GNP}/snp_vs_nonsnp/gea_classgain_fdr.csv", index=False)
print("built gain tables")
"""

md_t1 = r"""### Table 1 — Bonferroni (0.05) significant GEA peaks, per model (bio1, gen9)"""
code_t1 = "gea_bonf_df\n"
md_t2 = r"""### Table 2 — BH-FDR (q<0.05) significant GEA peaks, per model (bio1, gen9)"""
code_t2 = "gea_fdr_df\n"

md_genes = r"""## Genes under non-SNP-only / SV-only GEA peaks

For every clq0.9 block a class (non-SNP or SV) flags significant that the SNP scan does
not — in ANY model — map the block interval (+/-2 kb flank) to TAIR10 genes and record
which models flagged it. `klass` = layer, `overlap` = in_block vs flank2kb, `n_models` /
`models` = where the layer fired but SNPs didn't. Two thresholds (Bonferroni, FDR).
Subject to the per-locus null; hypothesis-generating."""

code_genes = r"""
FLANK = 2000
bmap = {r.block: (r.chrom, int(r.start_pos), int(r.end_pos)) for r in blocks_clq09.load_blocks().itertuples()}
genes = lib.load_genes()

# optional descriptions cache (offline-safe; no live Ensembl call)
DESC = "analysis/grenenet_gea/varexp/gene_descriptions.csv"
desc = {}
if os.path.exists(DESC):
    _d = pd.read_csv(DESC).fillna("")
    desc = {r.gene: (getattr(r, "symbol", ""), getattr(r, "function", "")) for r in _d.itertuples()}

def gene_table(thr):
    hits = {}
    for cls in ("nonsnp", "sv"):
        blk2models = {}
        for m in MODELS:
            only = S[(m, cls)][thr] - S[(m, "snp")][thr]
            for b in only:
                blk2models.setdefault(b, set()).add(m)
        for b, ms in blk2models.items():
            if b not in bmap:
                continue
            c, s0, e0 = bmap[b]
            inb = genes[(genes.chrom == c) & (genes.start <= e0) & (genes.end >= s0)]
            flk = genes[(genes.chrom == c) & (genes.start <= e0 + FLANK) & (genes.end >= s0 - FLANK)]
            flk = flk[~flk.gene.isin(inb.gene)]
            for gid, nm, otype in ([(g, n, "in_block") for g, n in zip(inb.gene, inb.name.fillna(""))]
                                   + [(g, n, "flank2kb") for g, n in zip(flk.gene, flk.name.fillna(""))]):
                d = hits.setdefault((gid, cls), dict(gene=gid, name=nm, klass=cls, chrom=c,
                                                     overlap=otype, blocks=set(), models=set()))
                d["blocks"].add(b); d["models"].update(ms)
                if otype == "in_block":
                    d["overlap"] = "in_block"
    rows = [dict(gene=d["gene"], name=d["name"],
                 symbol=desc.get(d["gene"], ("", ""))[0], function=desc.get(d["gene"], ("", ""))[1],
                 klass=d["klass"], chrom=d["chrom"], overlap=d["overlap"],
                 n_blocks=len(d["blocks"]), n_models=len(d["models"]),
                 models=";".join(sorted(d["models"]))) for d in hits.values()]
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["klass", "n_models", "overlap", "gene"],
                            ascending=[True, False, True, True]).reset_index(drop=True)
    return df

genes_bonf_df = gene_table(0)
genes_fdr_df = gene_table(1)
genes_bonf_df.to_csv(f"{GNP}/snp_vs_nonsnp/gea_gained_genes_bonf.csv", index=False)
genes_fdr_df.to_csv(f"{GNP}/snp_vs_nonsnp/gea_gained_genes_fdr.csv", index=False)
print(f"gained genes — Bonferroni: {len(genes_bonf_df)}   FDR: {len(genes_fdr_df)}")
"""

md_gb = r"""### Genes — Bonferroni (non-SNP-only / SV-only GEA peaks)"""
code_gb = "genes_bonf_df\n"
md_gf = r"""### Genes — BH-FDR (non-SNP-only / SV-only GEA peaks)"""
code_gf = "genes_fdr_df\n"

md_bottom = r"""## Bottom line
- `n_nonsnp_not_snp` / `n_sv_not_snp` = the blocks the non-SNP / SV layer adds beyond SNPs.
  In LFMM (the calibrated model) these are the defensible gains; binomial/Kendall gains are
  inflated by pseudoreplication and should be read as an upper bound.
- The gene tables list what those gained blocks contain — hypothesis-generating, per-locus
  null applies. This is the bio1 view; the per-axis extension (bio5/12/13/16/19/pc1) is the
  same tables with one row per (axis, model)."""

nb = new_notebook(cells=[
    new_markdown_cell(md_title), new_code_cell(code_load),
    new_markdown_cell(md_tables), new_code_cell(code_tables),
    new_markdown_cell(md_t1), new_code_cell(code_t1),
    new_markdown_cell(md_t2), new_code_cell(code_t2),
    new_markdown_cell(md_genes), new_code_cell(code_genes),
    new_markdown_cell(md_gb), new_code_cell(code_gb),
    new_markdown_cell(md_gf), new_code_cell(code_gf),
    new_markdown_cell(md_bottom),
])
ep = ExecutePreprocessor(timeout=1200, kernel_name="python3", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": HERE}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
