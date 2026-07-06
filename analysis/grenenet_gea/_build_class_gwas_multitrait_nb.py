#!/usr/bin/env python
"""Build+execute the MULTI-TRAIT (all-sites Bolormaa meta) class-split GWAS notebook
(user 2026-07-03). Focus: SNP-vs-nonSNP significant-BLOCK overlap (clq0.9 blocks, matching
phase1_replication/_build_manhattan_clq90_nb.py's overlap_df convention) for JOINT/GLOBAL, and
CLIMATE re-derived for all 19 bioclim variables + PC1 (cheap: reuses the already-saved per-site
Z matrices, no GWAS rerun). No Manhattans -- table-only, per user request. Companion:
_build_class_gwas_persite_nb.py (the per-garden view). Runs in `basic` env.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_gea/notebooks/class_gwas_multitrait.ipynb"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

md_title = r"""# Multi-trait (all-sites) founder GWAS — SNP vs non-SNP peak overlap

**This notebook is the MULTI-TRAIT view**: each marker's 31 per-site z-scores (one LOCO-EMMAX scan
per garden) are combined with a Bolormaa meta-analysis into three genome-wide contrasts —

- **JOINT** (31 df, χ²) — selected at *any* site
- **GLOBAL** (1 df) — generalist: same-direction effect across all 31 sites
- **CLIMATE** (1 df) — local adaptation: effect tracks a site climate variable

run on SNP-only markers (1.75M) vs non-SNP markers pooled (indel+SV, 536k) -- same trait, same
shared LOCO kinship correction, only the tested marker class differs.

**Table, not Manhattans**: for each contrast, markers are collapsed to their **clq0.9 LD block**
(`lib.assign_clq_blocks`, the same r²≥0.9 partition used in the phase-1-replication SNP-vs-nonSNP
overlap table), one lead marker per block, BH q<0.05 across blocks -- then we ask how much the
SNP-significant and non-SNP-significant block sets **overlap** (Jaccard etc.), exactly like the
kendall/lfmm/binomial table in `phase1_replication/_build_manhattan_clq90_nb.py` but with
JOINT/GLOBAL/CLIMATE as the rows instead of GEA methods.

**CLIMATE, all bioclim variables**: bio1 (mean annual temperature) was the only climate axis
tested before. Re-derived here for **all 19 bioclim variables + PC1 of all 19** -- cheaply, by
reusing the already-computed per-site Z matrices (the Bolormaa CLIMATE contrast is just a
different linear combination of the same 31 per-site z-scores; JOINT/GLOBAL don't depend on
climate and are unchanged). For the per-garden (31 individual scans) view, see
`class_gwas_persite.ipynb`. Full design in `VAREXP_SELECTION_HANDOFF.md`."""

code_load = r"""
import os, sys, json
import numpy as np, pandas as pd
from scipy import stats
os.chdir("/global/scratch/users/tbellg/kmate")
sys.path.insert(0, "analysis/grenenet_gea")
import lib
OUT = "results/grenenet_gea/varexp"
CLASSES = ["snp", "nonsnp"]

def load_class(name):
    d = np.load(f"{OUT}/class_gwas_{name}.npz", allow_pickle=True)
    return dict(chrom=d["chrom"], pos=d["pos"], Z=d["Z"], sites=d["sites"], bio1=d["bio1"],
               p_joint=d["p_joint"], p_global=d["p_global"])

data = {c: load_class(c) for c in CLASSES}
sites = data["snp"]["sites"]
print(f"{len(sites)} sites | SNP M={len(data['snp']['pos']):,}  non-SNP M={len(data['nonsnp']['pos']):,}")
"""

code_helpers = r"""
def block_ids(chrom_lower, pos):
    chrom_cap = np.array([c.replace("chr", "Chr") for c in chrom_lower])
    return lib.assign_clq_blocks(chrom_cap, pos, r2=0.9)

def collapse_lead(chrom, pos, p):
    '''Collapse markers to one LEAD (min-p) marker per clq0.9 block. Markers landing in an
    inter-block gap ('') are dropped. collapse_to_blocks selects by MAX ABS(stat), so pass
    -log10(p) (always >=0, bigger = more significant) rather than p or -p directly -- abs(-p)
    is just p again and picks the wrong end.'''
    blk = block_ids(chrom, pos)
    keep = blk != ""
    nlp = -np.log10(np.clip(p[keep], 1e-300, 1))
    df = pd.DataFrame({"block": blk[keep], "nlp": nlp})
    lead = lib.collapse_to_blocks(df, stat="nlp", block_col="block")
    lead["p"] = 10.0 ** (-lead["nlp"])
    return lead

def sig_sets(lead, q=0.05):
    '''Block-level significant sets at BOTH Bonferroni and BH-FDR thresholds, from one lead-marker
    collapse (avoids re-collapsing per threshold).'''
    bonf = 0.05 / len(lead)
    bonf_sig = set(lead.loc[lead["p"] < bonf, "block"])
    qv = lib.bh(lead["p"].to_numpy())
    fdr_sig = set(lead.loc[qv < q, "block"])
    return bonf_sig, fdr_sig

def overlap_row(name, snp_bonf, snp_fdr, nonsnp_bonf, nonsnp_fdr):
    '''Overlap (jaccard/pct) is computed at the FDR threshold (the more inclusive, standard
    comparison set); Bonferroni is reported as a stricter count alongside it.'''
    shared = snp_fdr & nonsnp_fdr; union = snp_fdr | nonsnp_fdr
    return dict(contrast=name, n_snp_bonf=len(snp_bonf), n_snp_fdr=len(snp_fdr),
               n_nonsnp_bonf=len(nonsnp_bonf), n_nonsnp_fdr=len(nonsnp_fdr),
               n_shared_fdr=len(shared),
               jaccard_fdr=round(len(shared) / len(union), 3) if union else np.nan,
               pct_snp_also_nonsnp=round(len(shared) / len(snp_fdr), 3) if snp_fdr else np.nan,
               pct_nonsnp_also_snp=round(len(shared) / len(nonsnp_fdr), 3) if nonsnp_fdr else np.nan)

# block ids computed once per class (shared across all contrasts below)
blk_cache = {c: block_ids(data[c]["chrom"], data[c]["pos"]) for c in CLASSES}
"""

md_climate_setup = """## CLIMATE, re-derived for all 19 bioclim variables + PC1

Cheap: the Bolormaa CLIMATE contrast is `z_clim = (Z @ Cinv @ c) / sqrt(c @ Cinv @ c)` where `c`
is the (GLOBAL-orthogonalized) standardized site climate vector -- swapping in a different bioclim
variable only changes `c`, not the underlying per-site Z or the JOINT/GLOBAL results. No GWAS
rerun needed."""

code_climate_setup = r"""
clim_all = lib.load_climate().reindex(sites)[lib.BIO_COLS]
assert not clim_all.isna().any().any(), "missing bioclim value for a GWAS site"
Bz = (clim_all - clim_all.mean()) / clim_all.std()
Uu, Ss, Vt = np.linalg.svd(Bz.to_numpy(), full_matrices=False)
pc1 = Uu[:, 0] * Ss[0]
climate_axes = {f"bio{i}": clim_all[f"bio{i}"].to_numpy(float) for i in range(1, 20)}
climate_axes["PC1_allbio"] = pc1
print(f"{len(climate_axes)} climate axes: {list(climate_axes)[:3]} ... PC1_allbio "
      f"(explains {Ss[0]**2/np.sum(Ss**2):.1%} of standardized bioclim variance)")

# per-class corrcoef(Z) computed once, reused across all 20 climate axes
Cmat = {c: np.corrcoef(data[c]["Z"].T) for c in CLASSES}

def climate_p_cached(cls, bvec):
    Z = data[cls]["Z"]; C = Cmat[cls]; S = Z.shape[1]
    Cinv = np.linalg.pinv(C); one = np.ones(S)
    dg = float(one @ Cinv @ one)
    c0 = (bvec - bvec.mean()) / bvec.std()
    c = c0 - (float(one @ Cinv @ c0) / dg) * one
    dc = float(c @ Cinv @ c)
    z_clim = (Z @ Cinv @ c) / np.sqrt(dc)
    return 2 * stats.norm.sf(np.abs(z_clim))
"""

md_table = """## The overlap table

Rows: JOINT, GLOBAL, and CLIMATE for each of the 19 bioclim variables + PC1 of all 19 (22 rows
total). Columns match the phase-1-replication kendall/lfmm/binomial overlap table, plus both
significance thresholds shown side by side: `n_*_bonf` (Bonferroni 0.05, strict) and `n_*_fdr`
(BH q<0.05, the threshold the overlap/jaccard/pct columns are computed on)."""

code_table = r"""
rows = []
for name, key in [("JOINT", "p_joint"), ("GLOBAL", "p_global")]:
    lead = {c: collapse_lead(data[c]["chrom"], data[c]["pos"], data[c][key]) for c in CLASSES}
    sig = {c: sig_sets(lead[c]) for c in CLASSES}
    rows.append(overlap_row(name, *sig["snp"], *sig["nonsnp"]))
for axname, bvec in climate_axes.items():
    lead = {c: collapse_lead(data[c]["chrom"], data[c]["pos"], climate_p_cached(c, bvec)) for c in CLASSES}
    sig = {c: sig_sets(lead[c]) for c in CLASSES}
    rows.append(overlap_row(f"CLIMATE_{axname}", *sig["snp"], *sig["nonsnp"]))
overlap_df = pd.DataFrame(rows)
overlap_df
"""

md_bottom = """## Bottom line

- **JOINT/GLOBAL**: read off the `n_*_bonf` / `n_*_fdr` / `jaccard_fdr` / `pct_*` columns
  directly -- the same overlap-table logic as the phase-1 kendall/lfmm/binomial comparison, now
  applied to the SNP-vs-nonSNP class split, with both significance thresholds shown.
- **CLIMATE across all 19 bioclim variables + PC1**: if the climate-gradient null established
  earlier (bio1) generalizes, expect **0 (or near-0) significant blocks in most/all CLIMATE rows,
  for both classes** -- a systematic sweep across the full bioclim panel (not just temperature)
  before concluding local adaptation is absent for this system."""

nb = new_notebook(cells=[
    new_markdown_cell(md_title),
    new_code_cell(code_load),
    new_code_cell(code_helpers),
    new_markdown_cell(md_climate_setup),
    new_code_cell(code_climate_setup),
    new_markdown_cell(md_table),
    new_code_cell(code_table),
    new_markdown_cell(md_bottom),
])
ep = ExecutePreprocessor(timeout=1200, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
