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
OUT = f"{ROOT}/analysis/grenenet_selection/notebooks/class_gwas_multitrait.ipynb"
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
sys.path.insert(0, "analysis/grenenet_selection")
import lib
OUT = "analysis/grenenet_selection/r3_persite_gwas/results/varexp"
CLASSES = ["snp", "nonsnp", "sv"]

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

md_overlap = """## Peak overlap — two tables (Bonferroni & FDR)

Peaks = **clq0.9-block leads significant** at each threshold (one lead marker per LD block, so
markers in LD aren't double-counted). We build the **same table twice** — once with peaks called at
**Bonferroni** 0.05 (strict) and once at **BH-FDR** q<0.05 (inclusive) — for all three classes.
Rows: JOINT, GLOBAL, and CLIMATE for each of the 19 bioclim variables + PC1 of all 19.

The headline is **what the non-SNP / SV layer adds that SNPs miss** — not the shared blocks. So the
key columns are `n_nonsnp_not_snp` / `n_sv_not_snp`: blocks that class flags as significant but the
SNP scan does **not** (blocks live on a shared clq0.9 partition, so the set difference is exact).
`pct_*_not_snp` is that count as a fraction of the class's own hits. `n_snp`/`n_nonsnp`/`n_sv` give
the raw per-class hit totals for context. A block counts as "not in snp" whether SNPs there are
merely non-significant *or* absent entirely — both are signal SNPs don't deliver."""

code_overlap_setup = r"""
def lead_from_p(c, p):
    '''Collapse to one lead (min-p) marker per clq0.9 block, reusing the cached block ids.'''
    b = blk_cache[c]; keep = b != ""
    nlp = -np.log10(np.clip(p[keep], 1e-300, 1))
    lead = lib.collapse_to_blocks(pd.DataFrame({"block": b[keep], "nlp": nlp}),
                                  stat="nlp", block_col="block")
    lead["p"] = 10.0 ** (-lead["nlp"])
    return lead

def overlap_row(name, sets):
    # The gain from the non-SNP / SV layer = blocks IT calls significant that the SNP scan does NOT
    # (blocks are on a shared clq0.9 partition, so the set difference is well-defined).
    row = {"contrast": name, "n_snp": len(sets["snp"]),
           "n_nonsnp": len(sets["nonsnp"]), "n_sv": len(sets["sv"])}
    for b in ("nonsnp", "sv"):
        new = sets[b] - sets["snp"]          # significant in b, NOT in snp = the new information
        row[f"n_{b}_not_snp"] = len(new)
        row[f"pct_{b}_not_snp"] = round(len(new) / len(sets[b]), 3) if sets[b] else np.nan
    return row

# per (contrast, class): compute the lead collapse once, keep BOTH (bonf, fdr) significant sets
contrasts = [("JOINT", lambda c: data[c]["p_joint"]),
             ("GLOBAL", lambda c: data[c]["p_global"])]
for axname, bvec in climate_axes.items():
    contrasts.append((f"CLIMATE_{axname}", lambda c, bv=bvec: climate_p_cached(c, bv)))

sig_by_contrast = {name: {c: sig_sets(lead_from_p(c, pf(c))) for c in CLASSES}
                   for name, pf in contrasts}

def overlap_table(thr_idx):   # 0 = Bonferroni, 1 = FDR
    return pd.DataFrame([overlap_row(name, {c: sig_by_contrast[name][c][thr_idx] for c in CLASSES})
                         for name, _ in contrasts])

overlap_bonf_df = overlap_table(0)
overlap_fdr_df = overlap_table(1)
print("built Bonferroni + FDR overlap tables:", overlap_bonf_df.shape, overlap_fdr_df.shape)
"""

md_bonf = """### Table 1 — Bonferroni (0.05) significant peaks"""
code_bonf = "overlap_bonf_df\n"

md_fdr = """### Table 2 — BH-FDR (q<0.05) significant peaks"""
code_fdr = "overlap_fdr_df\n"

md_genes = """## Genes under non-SNP-only / SV-only peaks — anything interesting?

For every block that is significant in **non-SNP** or in **SV** but **not** in SNP (the gain columns
above), map the block's clq0.9 interval (±2 kb promoter flank) to TAIR10 genes, and record **which
contrast(s)** flagged it — i.e. JOINT, GLOBAL, or which bioclim variable (bio1–19 / PC1). Two tables:
peaks called at **Bonferroni** and at **FDR**. `klass` = which layer (nonsnp/sv), `overlap` =
in_block vs flank2kb, `contrasts` = the bio variable(s) / omnibus it was significant for. This is a
hypothesis-generating browse (the per-locus climate null still applies), sorted so genes recurring
across the most contrasts surface first."""

code_genes_setup = r"""
FLANK = 2000
def block_interval_map():
    m = {}
    for ci in range(1, 6):
        c = f"Chr{ci}"
        g = pd.read_csv(f"{lib.CLQ_BLOCKS_DIR}/chr{ci}_clq0.9_blocks_clq0.9.tsv",
                        sep="\t").sort_values("start_pos").reset_index(drop=True)
        for idx, r in g.iterrows():
            m[f"{c}_{idx}"] = (c, int(r.start_pos), int(r.end_pos))
    return m
bmap = block_interval_map()
genes = lib.load_genes()

def gene_table(thr_idx):
    hits = {}
    for b_cls in ("nonsnp", "sv"):
        blk2contrasts = {}
        for name, _ in contrasts:
            only = sig_by_contrast[name][b_cls][thr_idx] - sig_by_contrast[name]["snp"][thr_idx]
            for b in only:
                blk2contrasts.setdefault(b, set()).add(name.replace("CLIMATE_", ""))
        for b, cs in blk2contrasts.items():
            if b not in bmap:
                continue
            c, s, e = bmap[b]
            inb = genes[(genes.chrom == c) & (genes.start <= e) & (genes.end >= s)]
            flk = genes[(genes.chrom == c) & (genes.start <= e + FLANK) & (genes.end >= s - FLANK)]
            flk = flk[~flk.gene.isin(inb.gene)]
            for gid, nm, otype in ([(g, n, "in_block") for g, n in zip(inb.gene, inb.name.fillna(""))]
                                   + [(g, n, "flank2kb") for g, n in zip(flk.gene, flk.name.fillna(""))]):
                d = hits.setdefault((gid, b_cls), dict(gene=gid, name=nm, klass=b_cls, chrom=c,
                                                       overlap=otype, blocks=set(), contrasts=set()))
                d["blocks"].add(b); d["contrasts"].update(cs)
                if otype == "in_block":
                    d["overlap"] = "in_block"
    rows = [dict(gene=d["gene"], name=d["name"], klass=d["klass"], chrom=d["chrom"], overlap=d["overlap"],
                 n_blocks=len(d["blocks"]), n_contrasts=len(d["contrasts"]),
                 contrasts=";".join(sorted(d["contrasts"]))) for d in hits.values()]
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["klass", "n_contrasts", "overlap", "gene"],
                            ascending=[True, False, True, True]).reset_index(drop=True)
    return df

genes_bonf_df = gene_table(0)
genes_fdr_df = gene_table(1)

# ---- Ensembl Plants gene descriptions (cached to gene_descriptions.csv; offline -> AT-IDs only) ----
import urllib.request as _u
DESC_CSV = f"{OUT}/gene_descriptions.csv"
_cache = {}
if os.path.exists(DESC_CSV):
    _c = pd.read_csv(DESC_CSV).fillna("")
    _cache = {r.gene: (r.symbol, r.function) for r in _c.itertuples()}
_need = sorted(set(genes_bonf_df.gene) | set(genes_fdr_df.gene))
_missing = [g for g in _need if g not in _cache]

def _ens_batch(ids, chunk=45):
    got = {}
    for i in range(0, len(ids), chunk):
        sub = ids[i:i + chunk]
        for att in range(4):
            try:
                req = _u.Request("https://rest.ensembl.org/lookup/id",
                                 data=json.dumps({"ids": sub}).encode(),
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
                for g, v in json.load(_u.urlopen(req, timeout=45)).items():
                    v = v or {}
                    got[g] = (v.get("display_name", "") or "",
                              (v.get("description", "") or "").split(" [Source")[0])
                break
            except Exception as e:
                import time as _t; _t.sleep(2 * (att + 1))
                if att == 3:
                    print(f"[warn] ensembl chunk {i} failed after retries: {type(e).__name__}")
    return got

if _missing:
    new = _ens_batch(_missing)
    _cache.update(new)
    if _cache:
        pd.DataFrame([{"gene": g, "symbol": s, "function": f} for g, (s, f) in sorted(_cache.items())]
                     ).to_csv(DESC_CSV, index=False)
    print(f"Ensembl: described {len(new)}/{len(_missing)} new genes "
          f"({'offline -> AT-IDs only' if not new else f'cache -> {len(_cache)} total'})")
else:
    print(f"Ensembl: all {len(_need)} genes from cache")

def _add_desc(df):
    df = df.copy()
    df.insert(2, "symbol", df.gene.map(lambda g: _cache.get(g, ("", ""))[0]))
    df.insert(3, "function", df.gene.map(lambda g: _cache.get(g, ("", ""))[1]))
    return df
genes_bonf_df = _add_desc(genes_bonf_df)
genes_fdr_df = _add_desc(genes_fdr_df)

print(f"non-SNP/SV-only genes — Bonferroni: {len(genes_bonf_df)}   FDR: {len(genes_fdr_df)}")
for lbl, df in [("Bonferroni", genes_bonf_df), ("FDR", genes_fdr_df)]:
    if len(df):
        print(f"  {lbl}: " + ", ".join(f"{k}={int(v)}" for k, v in df.klass.value_counts().items())
              + f"; in_block={(df.overlap=='in_block').sum()}")
"""

md_genes_bonf = """### Genes — Bonferroni (non-SNP-only / SV-only)"""
code_genes_bonf = "genes_bonf_df\n"
md_genes_fdr = """### Genes — FDR (non-SNP-only / SV-only)"""
code_genes_fdr = "genes_fdr_df\n"

md_bottom = """## Bottom line

- **Two tables, two thresholds**: `overlap_bonf_df` (strict) and `overlap_fdr_df` (inclusive). The
  quantity that matters is `n_nonsnp_not_snp` / `n_sv_not_snp` — blocks the non-SNP or SV scan finds
  that the SNP scan does **not**. That is the information gained by adding the layer; a small number
  (relative to `n_snp`) means the layer is mostly re-finding SNP peaks.
- **JOINT** is where hits concentrate; **GLOBAL** and every **CLIMATE** row (all 19 bioclim + PC1)
  should be ~empty for all three classes if the climate-gradient null holds — a full-panel sweep,
  not just temperature.
- Peaks are *significant* clq0.9 blocks (block-level BH over ~58k blocks, so FDR counts are much
  larger than the marker-level counts in `class_gwas_summary.json`) — peak-specific, not the
  whole-genome window concordance that was diluted by ~5.5k mostly-null windows.
- **`genes_bonf_df` / `genes_fdr_df`** list the genes under non-SNP-only / SV-only peaks and the bio
  variable(s) they were significant for — a candidate browse for kMate-unique signal, subject to the
  standing per-locus climate null."""

nb = new_notebook(cells=[
    new_markdown_cell(md_title),
    new_code_cell(code_load),
    new_code_cell(code_helpers),
    new_markdown_cell(md_climate_setup),
    new_code_cell(code_climate_setup),
    new_markdown_cell(md_overlap),
    new_code_cell(code_overlap_setup),
    new_markdown_cell(md_bonf),
    new_code_cell(code_bonf),
    new_markdown_cell(md_fdr),
    new_code_cell(code_fdr),
    new_markdown_cell(md_genes),
    new_code_cell(code_genes_setup),
    new_markdown_cell(md_genes_bonf),
    new_code_cell(code_genes_bonf),
    new_markdown_cell(md_genes_fdr),
    new_code_cell(code_genes_fdr),
    new_markdown_cell(md_bottom),
])
ep = ExecutePreprocessor(timeout=1200, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
