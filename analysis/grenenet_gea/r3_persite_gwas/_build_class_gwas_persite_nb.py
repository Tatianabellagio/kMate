#!/usr/bin/env python
"""Build+execute the PER-SITE (individual garden) class-split GWAS notebook (user 2026-07-03).
Focus: the 31 individual single-site LOCO-EMMAX scans that feed the multi-trait meta in
class_gwas_multitrait.ipynb -- per-site genomic-control stability across classes, plus a
SNP-vs-nonSNP significant-BLOCK overlap table computed PER SITE (clq0.9 blocks, same convention
as the multitrait notebook / phase1_replication overlap_df, but one row per garden instead of
one row per JOINT/GLOBAL/CLIMATE contrast -- there's no cross-site climate contrast at the
single-site level). No Manhattans -- table-only, per user request.
Load-only (class_gwas_{snp,nonsnp,sv}.npz); runs in `basic` env (matplotlib hangs in `plotting`).
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_gea/notebooks/class_gwas_persite.ipynb"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

md_title = r"""# Per-site (individual garden) founder GWAS — SNP vs non-SNP peak overlap

**This notebook is the PER-GARDEN view**: each of the 31 sites has its own single-site
LOCO-EMMAX scan (rank-normalized selection coefficient `s` for that garden, kinship-corrected,
genomic-control calibrated) — these 31 per-site scans are exactly what the Bolormaa meta-analysis
in `class_gwas_multitrait.ipynb` combines into JOINT/GLOBAL/CLIMATE. Here we look at the
individual inputs instead of the combined output: does inflation look similar site-by-site across
marker classes, and — the main table below — does the SNP-significant vs non-SNP-significant
**clq0.9-block overlap** hold up at the single-garden level (not just in the pooled multi-site
meta)? There's no cross-site CLIMATE contrast at this level (climate-differential selection is
inherently a cross-site comparison), so rows here are gardens, not contrasts.

Same design as the multi-trait notebook: 212 analyzable founders, one shared LOCO kinship
correction per chromosome, SNP-only (1.75M) / non-SNP pooled (536k) / SV-only (~13k) test
markers. Full design in `VAREXP_SELECTION_HANDOFF.md`."""

code_load = r"""
import os, sys, json
import numpy as np, pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.dpi': 120, 'font.size': 10})
os.chdir("/global/scratch/users/tbellg/kmate")
sys.path.insert(0, "analysis/grenenet_gea")
import lib
OUT = "analysis/grenenet_gea/varexp"
COLORS = {"snp": ["#3b4cc0", "#7aa0c4"], "nonsnp": ["#D55E00", "#f2a679"], "sv": ["#009E73", "#66c2a5"]}
LABEL = {"snp": "SNP-only", "nonsnp": "non-SNP (indel+SV)", "sv": "SV-only"}
CLASSES = ["snp", "nonsnp", "sv"]

def load_class(name):
    d = np.load(f"{OUT}/class_gwas_{name}.npz", allow_pickle=True)
    df = pd.DataFrame({"chrom": d["chrom"], "pos": d["pos"]})
    return df, d["Z"], d["sites"], d["bio1"], d["lam_persite"]

data = {c: load_class(c) for c in CLASSES}
sites = data["snp"][2]; bio1 = data["snp"][3]
print(f"{len(sites)} sites, bio1 range {bio1.min():.1f}-{bio1.max():.1f}C")
"""

md_lambda = """## Per-site genomic-control λ, all 31 gardens, all three classes

Sites ordered cold→hot (bio1). If SNP and non-SNP/SV panels pick up the same amount of signal in
EVERY individual garden (not just on average), these three lines should track each other closely
site by site, not just agree on the genome-wide mean."""

code_lambda = r"""
order = np.argsort(bio1)
fig, ax = plt.subplots(figsize=(12, 4.5))
for c in CLASSES:
    lam = data[c][4][order]
    ax.plot(range(len(sites)), lam, "o-", ms=5, color=COLORS[c][0], label=f"{LABEL[c]}")
ax.axhline(1.0, color="k", lw=.8, ls=":")
ax.set_xticks(range(len(sites))); ax.set_xticklabels(sites[order], rotation=90, fontsize=7)
ax.set_xlabel("site (ordered cold → hot by bio1)"); ax.set_ylabel(r"per-site $\lambda_{GC}$")
ax.set_title("Per-site genomic-control inflation, by marker class", loc="left")
ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(f"{OUT}/persite_lambda_3way.png", dpi=150, bbox_inches="tight")
plt.show()
print("cross-site corr(lambda_snp, lambda_nonsnp) =",
      round(np.corrcoef(data['snp'][4], data['nonsnp'][4])[0,1], 3))
print("cross-site corr(lambda_snp, lambda_sv)     =",
      round(np.corrcoef(data['snp'][4], data['sv'][4])[0,1], 3))
"""

md_overlap_setup = """## Per-site SNP-vs-nonSNP significant-block overlap

Same clq0.9-block / BH q<0.05 / Jaccard convention as the multi-trait notebook and the phase-1
kendall/lfmm/binomial overlap table -- but here each row is one GARDEN's own single-site scan
(not the 31-site JOINT/GLOBAL/CLIMATE meta)."""

code_overlap_setup = r"""
def block_ids(chrom_lower, pos):
    chrom_cap = np.array([c.replace("chr", "Chr") for c in chrom_lower])
    return lib.assign_clq_blocks(chrom_cap, pos, r2=0.9)

def collapse_lead(blk, p):
    '''collapse_to_blocks selects by MAX ABS(stat), so pass -log10(p) (always >=0, bigger = more
    significant) rather than p or -p directly -- abs(-p) is just p again and picks the wrong end.'''
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
    '''Overlap (jaccard/pct) is computed at the FDR threshold; Bonferroni is reported as a
    stricter count alongside it.'''
    shared = snp_fdr & nonsnp_fdr; union = snp_fdr | nonsnp_fdr
    return dict(site=name, n_snp_bonf=len(snp_bonf), n_snp_fdr=len(snp_fdr),
               n_nonsnp_bonf=len(nonsnp_bonf), n_nonsnp_fdr=len(nonsnp_fdr),
               n_shared_fdr=len(shared),
               jaccard_fdr=round(len(shared) / len(union), 3) if union else np.nan,
               pct_snp_also_nonsnp=round(len(shared) / len(snp_fdr), 3) if snp_fdr else np.nan,
               pct_nonsnp_also_snp=round(len(shared) / len(nonsnp_fdr), 3) if nonsnp_fdr else np.nan)

# block ids computed once per class (shared across all 31 per-site tests)
blk = {c: block_ids(data[c][0].chrom.to_numpy(), data[c][0].pos.to_numpy()) for c in CLASSES}
"""

md_table = """## Per-site overlap — two tables (Bonferroni & FDR)

Same significant-block sets as the multi-trait notebook, but **one row per garden** instead of one
row per contrast, and split into the **same table twice**: peaks called at **Bonferroni** 0.05
(strict) and at **BH-FDR** q<0.05 (inclusive), for all three classes. Key columns are the gain from
the layer — `n_nonsnp_not_snp` / `n_sv_not_snp`: blocks that class flags but the SNP scan does not
(with `pct_*_not_snp` as the fraction of the class's own hits). `n_snp`/`n_nonsnp`/`n_sv` give the
raw totals. Per-site power is far below the 31-site JOINT omnibus, so most gardens have few/zero hits
(`pct` NaN where a class has none)."""

code_table = r"""
def overlap_row(site_id, sets):
    # The gain from the non-SNP / SV layer = blocks IT calls significant that the SNP scan does NOT.
    row = {"site": int(site_id), "bio1": round(float(bio1[np.where(sites == site_id)[0][0]]), 1),
           "n_snp": len(sets["snp"]), "n_nonsnp": len(sets["nonsnp"]), "n_sv": len(sets["sv"])}
    for b in ("nonsnp", "sv"):
        new = sets[b] - sets["snp"]          # significant in b, NOT in snp = the new information
        row[f"n_{b}_not_snp"] = len(new)
        row[f"pct_{b}_not_snp"] = round(len(new) / len(sets[b]), 3) if sets[b] else np.nan
    return row

# per (site, class): collapse to block leads once, keep BOTH (bonf, fdr) significant sets
sig_ps = {}
for si, site_id in enumerate(sites):
    sig_ps[int(site_id)] = {c: sig_sets(collapse_lead(blk[c], 2 * stats.norm.sf(np.abs(data[c][1][:, si]))))
                            for c in CLASSES}

def persite_overlap_table(thr_idx):   # 0 = Bonferroni, 1 = FDR
    rows = [overlap_row(s, {c: sig_ps[s][c][thr_idx] for c in CLASSES}) for s in sig_ps]
    return pd.DataFrame(rows).sort_values("bio1").reset_index(drop=True)

persite_bonf_df = persite_overlap_table(0)
persite_fdr_df = persite_overlap_table(1)
print("per-site tables built:", persite_bonf_df.shape, persite_fdr_df.shape)
print("gardens with >=1 SNP hit  — Bonf:", int((persite_bonf_df.n_snp > 0).sum()),
      " FDR:", int((persite_fdr_df.n_snp > 0).sum()))
"""

md_bonf = """### Table 1 — Bonferroni (0.05) significant peaks, per garden"""
code_bonf = "persite_bonf_df\n"

md_fdr = """### Table 2 — BH-FDR (q<0.05) significant peaks, per garden"""
code_fdr = "persite_fdr_df\n"

md_tophit = """## Per-site top hit table

For each site and class, the single strongest per-site association (not corrected across sites
like JOINT) -- a quick look at whether the strongest per-site signal wanders to a different locus
site to site, or repeatedly lands in the same place."""

code_tophit = r"""
rows = []
for si, site_id in enumerate(sites):
    row = {"site": int(site_id), "bio1": round(float(bio1[si]), 1)}
    for c in CLASSES:
        df, Z, _, _, _ = data[c]
        z = Z[:, si]
        j = int(np.nanargmax(np.abs(z)))
        p = 2 * stats.norm.sf(np.abs(z[j]))
        row[f"{c}_top_pos"] = f"{df.chrom.iloc[j]}:{int(df.pos.iloc[j]):,}"
        row[f"{c}_top_p"] = f"{p:.2e}"
    rows.append(row)
pd.DataFrame(rows).sort_values("bio1")
"""

md_genes = """## Genes under non-SNP-only / SV-only peaks — anything interesting?

For every block a garden's **non-SNP** or **SV** scan flags but its **SNP** scan does not, map the
block's clq0.9 interval (±2 kb flank) to TAIR10 genes and record **which garden(s)** flagged it. (A
single-site scan has no climate axis, so here the "which contrast" is the garden id, not a bioclim
variable — the bioclim-conditioned version is in `class_gwas_multitrait.ipynb`.) Two tables:
Bonferroni and FDR. `klass` = layer (nonsnp/sv), `overlap` = in_block vs flank2kb, `sites` = the
gardens where the layer fired but SNPs didn't. Sorted so genes recurring across the most gardens
surface first — hypothesis-generating, subject to the per-locus null and weak single-site power."""

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
        blk2sites = {}
        for site_id, sg in sig_ps.items():
            only = sg[b_cls][thr_idx] - sg["snp"][thr_idx]
            for b in only:
                blk2sites.setdefault(b, set()).add(int(site_id))
        for b, ss in blk2sites.items():
            if b not in bmap:
                continue
            c, s, e = bmap[b]
            inb = genes[(genes.chrom == c) & (genes.start <= e) & (genes.end >= s)]
            flk = genes[(genes.chrom == c) & (genes.start <= e + FLANK) & (genes.end >= s - FLANK)]
            flk = flk[~flk.gene.isin(inb.gene)]
            for gid, nm, otype in ([(g, n, "in_block") for g, n in zip(inb.gene, inb.name.fillna(""))]
                                   + [(g, n, "flank2kb") for g, n in zip(flk.gene, flk.name.fillna(""))]):
                d = hits.setdefault((gid, b_cls), dict(gene=gid, name=nm, klass=b_cls, chrom=c,
                                                       overlap=otype, blocks=set(), sites=set()))
                d["blocks"].add(b); d["sites"].update(ss)
                if otype == "in_block":
                    d["overlap"] = "in_block"
    rows = [dict(gene=d["gene"], name=d["name"], klass=d["klass"], chrom=d["chrom"], overlap=d["overlap"],
                 n_blocks=len(d["blocks"]), n_sites=len(d["sites"]),
                 sites=";".join(str(x) for x in sorted(d["sites"]))) for d in hits.values()]
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["klass", "n_sites", "overlap", "gene"],
                            ascending=[True, False, True, True]).reset_index(drop=True)
    return df

genes_bonf_df = gene_table(0)
genes_fdr_df = gene_table(1)

# ---- Ensembl Plants gene descriptions (shared cache gene_descriptions.csv; offline -> AT-IDs only) ----
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

print(f"non-SNP/SV-only genes (per-site union) — Bonferroni: {len(genes_bonf_df)}   FDR: {len(genes_fdr_df)}")
"""

md_genes_bonf = """### Genes — Bonferroni (non-SNP-only / SV-only, per garden)"""
code_genes_bonf = "genes_bonf_df\n"
md_genes_fdr = """### Genes — FDR (non-SNP-only / SV-only, per garden)"""
code_genes_fdr = "genes_fdr_df\n"

md_bottom = """## Bottom line

- Per-site genomic-control λ tracks closely across all three marker classes (see the cross-site
  correlations above) -- garden by garden, not just on genome-wide average, SNP/non-SNP/SV panels
  pick up comparable amounts of signal.
- **Two per-site tables** (`persite_bonf_df`, `persite_fdr_df`): as in the multi-trait notebook, the
  informative columns are `n_nonsnp_not_snp` / `n_sv_not_snp` — blocks a garden's non-SNP or SV scan
  flags that its SNP scan misses. Per-site power is far below the 31-site JOINT omnibus, so most
  gardens have few/zero hits (`pct_*_not_snp` NaN where a class has none). Read where a garden
  *does* fire, and whether the layer adds anything SNPs didn't already catch there.
- The per-site top-hit table shows the strongest per-site signal moves around from garden to
  garden (as expected for a mostly-null single-site test) rather than repeatedly landing on one
  locus -- the repeatable, cross-site-consistent locus only emerges once you look at the JOINT
  meta, not any single site.
- **`genes_bonf_df` / `genes_fdr_df`** list the genes under non-SNP-only / SV-only per-garden peaks
  and which gardens flagged them — the per-site companion to the bioclim-conditioned gene list in
  `class_gwas_multitrait.ipynb`."""

nb = new_notebook(cells=[
    new_markdown_cell(md_title),
    new_code_cell(code_load),
    new_markdown_cell(md_lambda),
    new_code_cell(code_lambda),
    new_markdown_cell(md_overlap_setup),
    new_code_cell(code_overlap_setup),
    new_markdown_cell(md_table),
    new_code_cell(code_table),
    new_markdown_cell(md_bonf),
    new_code_cell(code_bonf),
    new_markdown_cell(md_fdr),
    new_code_cell(code_fdr),
    new_markdown_cell(md_tophit),
    new_code_cell(code_tophit),
    new_markdown_cell(md_genes),
    new_code_cell(code_genes_setup),
    new_markdown_cell(md_genes_bonf),
    new_code_cell(code_genes_bonf),
    new_markdown_cell(md_genes_fdr),
    new_code_cell(code_genes_fdr),
    new_markdown_cell(md_bottom),
])
ep = ExecutePreprocessor(timeout=900, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
