#!/usr/bin/env python
"""Build + execute notebooks/nonsnp_specific_genes.ipynb (run in `basic` env).

Presents the nonSNP-specific climate-hit gene list (multiaxis/nonsnp_specific_genes.csv):
blocks BH-sig in non-SNP but NOT SNP across the 20 axes (bio1..bio19 + pc1) x 3 models,
with each gene's symbol + full functional DESCRIPTION, provenance (which axes/models),
and a functional-category tagging (heat / drought-ABA / cold / flowering-circadian).

  BPY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
  $BPY analysis/grenenet_gea/phase1_replication/multiaxis/_build_nonsnp_genes_nb.py
"""
import os
import nbformat as nbf
from nbconvert.preprocessors import ExecutePreprocessor

HERE = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication"
OUTDIR = f"{HERE}/notebooks"
OUT = f"{OUTDIR}/nonsnp_specific_genes.ipynb"
os.makedirs(OUTDIR, exist_ok=True)

nb = nbf.v4.new_notebook()
C = []

C.append(nbf.v4.new_markdown_cell(r"""# nonSNP-specific climate-hit genes (clq0.9, all 20 axes)

Genes on clq0.9 blocks that are WZA **BH q<0.05 in non-SNP (SV+small-indel) but NOT
in SNP**, for at least one of the 20 climate axes (bio1..bio19 + PC1) x 3 models
(kendall / lfmm-no-gif / quasi-binomial). These are regions where structural /
indel variation carries a climate association that SNPs miss.

Source: `multiaxis/nonsnp_specific_genes.csv` (block -> clq0.9 span -> overlapping
TAIR10 genes -> Ensembl Plants symbol + description). Provenance columns:
- `n_flags` = # of (axis,model) combos flagging the block nonSNP-specific
- `axes_nonsnp_specific`, `models` = which axes / models
- `never_snp_hit` = block is NEVER a SNP BH-hit on any axis/model (strongest nonSNP-only evidence)

**Caveat:** the pooled set includes the (uncalibrated, inflated) Kendall model; the
`models` column lets you restrict to the calibrated models (lfmm / binomial). The
`never_snp_hit=True` + high-`n_flags` rows are robust regardless."""))

C.append(nbf.v4.new_code_cell(r"""import pandas as pd, numpy as np
pd.set_option("display.max_rows", 600)
pd.set_option("display.max_colwidth", 90)
pd.set_option("display.width", 200)

CSV = "/global/scratch/users/tbellg/kmate/results/grenenet_gea/phase1_replication/multiaxis/nonsnp_specific_genes.csv"
d = pd.read_csv(CSV)
# one row per gene (a gene can appear via >1 block; keep the block with most flags)
g = (d[d.gene != ""].sort_values("n_flags", ascending=False)
     .drop_duplicates("gene").reset_index(drop=True))
g["description"] = g["description"].fillna(""); g["symbol"] = g["symbol"].fillna("")
print(f"{d.block.nunique()} nonSNP-specific blocks | {len(g)} unique genes "
      f"({(g.symbol!='').sum()} with symbols) | "
      f"{d[d.never_snp_hit==True].block.nunique()} blocks never a SNP hit")
g[["block","symbol","gene","description","n_flags","models","never_snp_hit"]].head(15)"""))

C.append(nbf.v4.new_markdown_cell(r"""## Functional-category tagging

Tag each gene by keyword in its symbol+description, plus a curated exact-symbol
list for well-known Arabidopsis genes whose description is terse."""))

C.append(nbf.v4.new_code_cell(r"""CATS = {
 "heat":       r"heat|thermo|high temperature|chaperone|heat shock|hsp|hsf|dnaj",
 "drought_ABA":r"drought|abscisic|\baba\b|dehydrat|desicc|osmotic|water stress|stomat|dehydrin|late embryogenesis|proline|salt",
 "cold":       r"\bcold\b|freezing|chilling|\bcbf\b|frost|cor15",
 "flower_circ":r"flower|floral|photoperiod|vernaliz|infloresc|florigen|circadian|clock|rhythm|oscillat",
}
CURATED = {  # exact symbol -> category (catches terse descriptions)
 "HSFA2":"heat","HSBP":"heat","HSP17.4":"heat","MBF1C":"heat",
 "RAS1":"drought_ABA","DREB2A":"drought_ABA","RD29A":"drought_ABA","NCED3":"drought_ABA",
 "COR15A":"cold","VRN2":"flower_circ","FT":"flower_circ","FLC":"flower_circ",
 "GI":"flower_circ","CCA1":"flower_circ","TOC1":"flower_circ","PIF4":"flower_circ",
}
txt = (g.symbol + " " + g.description).str.lower()
g["category"] = ""
for cat, pat in CATS.items():
    g.loc[(g.category == "") & txt.str.contains(pat, regex=True, na=False), "category"] = cat
for sym, cat in CURATED.items():
    g.loc[g.symbol == sym, "category"] = cat
print("category counts:\n", g["category"].replace("", "other").value_counts().to_string())"""))

C.append(nbf.v4.new_markdown_cell(r"""## Candidate genes by functional category (stress / phenology)

The scientifically interesting subset. HEAT dominates; flowering/circadian is nearly absent."""))

C.append(nbf.v4.new_code_cell(r"""for cat in ["heat", "drought_ABA", "cold", "flower_circ"]:
    sub = g[g.category == cat].sort_values("n_flags", ascending=False)
    print(f"\n===== {cat}: {len(sub)} genes =====")
    if len(sub):
        print(sub[["block","symbol","gene","description","n_flags","axes_nonsnp_specific","never_snp_hit"]]
              .to_string(index=False))"""))

C.append(nbf.v4.new_markdown_cell(r"""## Full nonSNP-specific gene list (with descriptions)

All unique genes, sorted by how many axis-model combos flag them. Written also to
`multiaxis/nonsnp_specific_genes_table.csv` (one row per gene)."""))

C.append(nbf.v4.new_code_cell(r"""full = g[["block","chrom","start","end","symbol","gene","description",
                          "category","n_flags","axes_nonsnp_specific","models","never_snp_hit"]].copy()
full.to_csv("/global/scratch/users/tbellg/kmate/results/grenenet_gea/phase1_replication/"
            "multiaxis/nonsnp_specific_genes_table.csv", index=False)
full"""))

nb["cells"] = C
nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3"}
print("executing…")
ep = ExecutePreprocessor(timeout=600, kernel_name="python3", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": HERE}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print("wrote", OUT)
