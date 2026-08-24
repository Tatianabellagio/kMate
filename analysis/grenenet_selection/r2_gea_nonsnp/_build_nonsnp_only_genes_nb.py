#!/usr/bin/env python
"""Build+execute the non-SNP-only candidate-genes notebook (user 2026-07-03): a readable annotated
table of the genes under clq0.9 blocks the non-SNP (indel+SV) GWAS flags but the SNP GWAS misses,
across all 53 contrasts (multitrait JOINT/GLOBAL + CLIMATE x [bio1-19, PC1] + 31 per-site), with a
+/-2kb promoter flank, Ensembl/UniProt descriptions, a Bonferroni-stricter subset, and a keyword
screen for flowering-time / cold / heat / circadian genes. Runs in `basic` env.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_selection/notebooks/nonsnp_only_genes.ipynb"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

md_title = r"""# Non-SNP-only candidate genes — what does kMate's indel+SV layer flag that SNPs miss?

For every contrast we ran (multitrait **JOINT**, **GLOBAL**, **CLIMATE** × [bio1–19 + PC1 of all
bioclim], and 31 **per-site** scans), we took the clq0.9 LD blocks BH-FDR-significant in the
**non-SNP (indel+SV)** GWAS but **not** in the SNP GWAS — the peaks a SNP-only GWAS would miss —
then mapped each block to overlapping **TAIR10 genes** (within the block span, plus a **±2 kb
promoter flank** to catch nearby regulatory targets) and described each gene via the **Ensembl
Plants** and **UniProt** REST APIs.

**Caveat up front (read before interpreting):** these are, by construction, the *fragile* end of
the signal — blocks one marker class calls significant and the other doesn't. The variance-
partition work showed non-SNP adds ~nothing genome-wide; this list is the scattered local
exceptions. Most are FDR-level (not Bonferroni), and the CLIMATE-contrast subset attaches to
bioclim hits that are likely FDR-tail noise (they mostly vanish under Bonferroni, don't replicate
across classes, and aren't corrected across the 20 collinear climate axes). Treat as a hypothesis-
generating candidate list, not confirmed kMate-unique adaptation loci."""

code_load = r"""
import os
import numpy as np, pandas as pd
os.chdir("/global/scratch/users/tbellg/kmate")
OUT = "analysis/grenenet_selection/varexp"
g = pd.read_csv(f"{OUT}/nonsnp_only_genes_described.csv").fillna("")
b = pd.read_csv(f"{OUT}/nonsnp_only_blocks.csv").fillna("")
g["desc"] = g["ensembl_description"].str.replace(r" \[Source:.*", "", regex=True)

def ctype(s):
    ks = set()
    for c in s.split(";"):
        if c.startswith("multitrait_CLIMATE"): ks.add("CLIMATE")
        elif c.startswith("multitrait_JOINT"): ks.add("JOINT")
        elif c.startswith("multitrait_GLOBAL"): ks.add("GLOBAL")
        elif c.startswith("persite"): ks.add("per-site")
    return ",".join(sorted(ks))
g["contrast_types"] = g["contrasts"].apply(ctype)
print(f"{len(g)} non-SNP-only candidate genes across {len(b)} blocks")
print(f"  in-block: {(g.overlap_type=='in_block').sum()}   +/-2kb flank-only: {(g.overlap_type=='flank2kb').sum()}")
print(f"  Bonferroni-subset genes: {g.bonferroni.sum()}")
print(f"  with an Ensembl description: {(g.ensembl_description!='').sum()};  with a UniProt function: {(g.uniprot_function!='').sum()}")
"""

md_bonf = """## Bonferroni-stricter subset — the most defensible non-SNP-only hits

Genes under blocks that are non-SNP-only at the **block-level Bonferroni** threshold (not just
FDR) in at least one contrast. This is the subset that survives the strict multiple-testing bar."""

code_bonf = r"""
bonf = g[g.bonferroni].copy().sort_values(["contrast_types", "gene"])
with pd.option_context("display.max_colwidth", 80, "display.width", 240):
    print(bonf[["gene", "symbol", "desc", "overlap_type", "contrast_types", "contrasts"]].to_string(index=False))
"""

md_theme = """## Themed screen — flowering time / cold / heat / circadian

Keyword match over each gene's symbol + Ensembl description + UniProt function text. This is
**high-precision / low-recall**: it catches genes whose annotation *explicitly* names the process,
but a gene with only a generic family name (e.g. "NAC domain protein") won't match even if it is
in fact stress-related. So absence here is not evidence of absence — see the manual notes below."""

code_theme = r"""
THEMES = {
    "flowering_time": ["flowering", "floral", "vernaliz", "photoperiod", "inflorescence",
                       "constans", "gigantea", "frigida", "agamous", "apetala", "flc", "soc1",
                       " ft ", "mads", "meristem identity"],
    "cold": ["cold", "freezing", "chilling", "cbf", "dreb", "cor15", "low temperature",
             "ice1", "dehydrin", "cold acclimation", "cold-regulated", "cold regulated"],
    "heat": ["heat shock", "heat stress", "high temperature", "thermotoler", "thermomorph",
             "hsp", "hsf", "chaperone"],
    "circadian": ["circadian", "clock", "cca1", "lhy", "toc1", "pseudo-response regulator",
                  "zeitlupe", "rhythm", " prr", "elf3", "elf4"],
}
text = (g["symbol"].str.lower() + " | " + g["desc"].str.lower() + " | " + g["uniprot_function"].str.lower())
for th, kws in THEMES.items():
    g[th] = text.apply(lambda t: any(k in t for k in kws))
g["themes"] = g.apply(lambda r: ",".join(th for th in THEMES if r[th]), axis=1)
hit = g[g["themes"] != ""]
print(f"{len(hit)} genes matched a theme keyword:")
with pd.option_context("display.max_colwidth", 90, "display.width", 250):
    print(hit[["gene", "symbol", "themes", "desc", "contrast_types"]].to_string(index=False))
"""

md_theme_fn = """### UniProt function text for the theme-matched genes (fuller context)"""

code_theme_fn = r"""
with pd.option_context("display.max_colwidth", 200, "display.width", 260):
    for _, r in hit.iterrows():
        fn = r["uniprot_function"] or "(no UniProt function annotation)"
        print(f"- {r['gene']} ({r['symbol'] or r['desc']}) [{r['themes']}]: {fn}\n")
"""

md_full = """## Full annotated table (all 126 genes)

Sorted Bonferroni-first, then by number of contrasts. `overlap_type` = in_block vs ±2kb flank."""

code_full = r"""
show = g.sort_values(["bonferroni", "n_contrasts", "gene"], ascending=[False, False, True])
with pd.option_context("display.max_rows", 200, "display.max_colwidth", 60, "display.width", 260):
    print(show[["gene", "symbol", "desc", "biotype", "overlap_type", "bonferroni",
                "n_contrasts", "contrast_types"]].to_string(index=False))
"""

md_bottom = """## Notes

- The **JOINT** (any-site selection) genes are the strongest-motivated subset (JOINT is the one
  contrast with a real, well-calibrated genome-wide signal); CLIMATE-contrast genes inherit the
  bioclim-null caveat above.
- Keyword screening is low-recall. A manual pass over the symbols (below the auto-screen) is worth
  doing for canonical stress/flowering/clock genes that carry only generic family annotations.
- Files: `analysis/grenenet_selection/varexp/nonsnp_only_{blocks,genes,genes_described}.csv`."""

nb = new_notebook(cells=[
    new_markdown_cell(md_title),
    new_code_cell(code_load),
    new_markdown_cell(md_bonf),
    new_code_cell(code_bonf),
    new_markdown_cell(md_theme),
    new_code_cell(code_theme),
    new_markdown_cell(md_theme_fn),
    new_code_cell(code_theme_fn),
    new_markdown_cell(md_full),
    new_code_cell(code_full),
    new_markdown_cell(md_bottom),
])
ep = ExecutePreprocessor(timeout=600, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
