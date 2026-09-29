#!/usr/bin/env python
"""Build SITE-COLLAPSED raw mirror Manhattans, LFMM only, all 20 climate axes — clq0.9
TILING partition. Site twin of `_build_raw_manhattan_lfmm_tile_nb.py`.

Difference from the pool notebook: the per-record LFMM p is computed on the SITE-MEAN Δp
(mean across the ~11 pools within each of the 31 sites), lfmm_ridge(K=1) + lfmm_test.
Collapsing pools -> sites removes the pseudoreplication that inflates the pool-level scan
(pool λ≈2.5 -> site λ≈1.7-1.8), and it acts as a filter: pool peaks driven by a few
correlated pools within a site collapse, while site-coherent climate peaks persist.
Inputs built by wza_investigation/deinflate/build_site_wza_in.py ->
results/multiaxis/wza_in_clq09_tile_site/lfmm_{cls}_gen9_{axis}.csv (pval = site RAW p).

Same "SNP up / class down" mirror-Manhattan + QQ style as the pool notebook.
Emits: raw_manhattan_snp_vs_sv_lfmm_site_tile.ipynb
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = "lfmm"


def fdr_md(cls, title_word):
    return f"""## FDR-relaxed gene list — {title_word}-only new peaks at BH **q < 0.05** (site-collapsed)

The mirror-Manhattans above call new peaks at the stringent per-record **Bonferroni** bar
(0.05/n). Here we relax to a **Benjamini–Hochberg FDR** of **q < 0.05**, computed on the same
site-level RAW p (per class, genome-wide within each axis). A "{cls}-only FDR new peak" = a
tiling block whose lead {cls} record clears q<0.05 while its lead SNP record does not — the FDR
analogue of the Bonferroni new-peak rule. Consolidated across all 20 axes, deduplicated per gene
(kept: strongest {cls} q), TAIR/UniProt-annotated, written to
`raw_manhattan_newpeak_genes_{cls}_site_tile_FDR.csv`.

**Caveat:** FDR is applied to the RAW (uncalibrated) p, so on the inflated axes (site λ up to
~3 on bio15/bio17) this list will inherit that inflation — read the per-axis λ in the QQ panels
above when weighing these genes."""


# self-contained FDR cell (redefines its own helpers so it can be appended without re-running
# the manhattan cells). __CLS__ is substituted per class to dodge f-string brace escaping.
_FDR_CODE = r'''
import os, sys, importlib.util as _ilu
import numpy as np, pandas as pd
sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), "..", "..")))
import lib
GEA = lib.GEA
WZAIN = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/wza_in_clq09_tile_site"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
CLS = "__CLS__"; QCUT = 0.05

def _spans(r2=0.9):
    tag = f"clq{r2}"; rows = []
    for ci in range(1, 6):
        ch = f"Chr{ci}"
        g = pd.read_csv(f"{lib.CLQ_BLOCKS_DIR}/chr{ci}_{tag}_blocks_{tag}.tsv", sep="\t").sort_values("start_pos").reset_index(drop=True)
        ends = g.end_pos.to_numpy(np.int64)
        for i in range(len(g)):
            lo = 1 if i == 0 else int(ends[i-1]) + 1
            hi = int(ends[i]) if i < len(g)-1 else 10**9
            rows.append((f"{ch}_{i}", ch, lo, hi))
    return pd.DataFrame(rows, columns=["block", "chrom", "start", "end"]).set_index("block")
SPANS = _spans(); GENES = lib.load_genes()
def genes_on(b):
    if b not in SPANS.index: return []
    s = SPANS.loc[b]; gc = GENES[(GENES.chrom == s.chrom) & (GENES.end >= s.start) & (GENES.start <= s.end)]
    return list(gc.gene)

def bh_q(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    ranked = (p[o] * n) / np.arange(1, n + 1)
    qo = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n); q[o] = np.clip(qo, 0, 1); return q

def load_q(cls, axis):
    f = f"{WZAIN}/lfmm_{cls}_gen9_{axis}.csv"
    if not os.path.exists(f): return None
    d = pd.read_csv(f); d = d[d["MAF"] > 0.05].copy()
    d["q"] = bh_q(d["pval"].to_numpy()); d["nlp"] = -np.log10(d["pval"].clip(lower=1e-300))
    return d

def leads_q(d):
    d = d[d["block"].notna() & (d["block"] != "")]
    return d.loc[d.groupby("block")["q"].idxmin()]

rows = []
for axis in AXES:
    s = load_q("snp", axis); c = load_q(CLS, axis)
    if s is None or c is None:
        print(f"missing {axis}"); continue
    ls_ = leads_q(s)[["block", "q"]].rename(columns={"q": "q_snp"})
    lc_ = leads_q(c)
    m = lc_.merge(ls_, on="block", how="left"); m["q_snp"] = m["q_snp"].fillna(1.0)
    nk = m[(m["q"] < QCUT) & (m["q_snp"] >= QCUT)].copy()
    nk["genes"] = nk["block"].map(lambda b: ";".join(genes_on(b))); nk["axis"] = axis
    rows.append(nk[["block", "chrom", "pos", "q", "nlp", "genes", "axis"]])
    print(f"{axis}: {len(nk)} {CLS}-only FDR(q<{QCUT}) new peaks (gene-bearing {(nk.genes != '').sum()})")
ALLF = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["block","chrom","pos","q","nlp","genes","axis"])

gr = []
for _, r in ALLF.iterrows():
    for g in (r.genes.split(";") if r.genes else []):
        gr.append(dict(gene=g, block=r.block, chrom=r.chrom, pos=r.pos, axis=r.axis, q=r.q, nlp=r.nlp))
GF = pd.DataFrame(gr)
if len(GF):
    recur = GF.groupby("gene")["axis"].nunique().rename("n_axes")
    best = GF.sort_values("q").drop_duplicates("gene").set_index("gene")
    candf = best.join(recur).sort_values(["n_axes", "q"], ascending=[False, True]).reset_index()
    _annp = os.path.abspath(os.path.join(os.getcwd(), "..", "annotate_genes_tair_uniprot.py"))
    _sp = _ilu.spec_from_file_location("ann_tu", _annp); _ann = _ilu.module_from_spec(_sp); _sp.loader.exec_module(_ann)
    _A = _ann.annotate(sorted(candf["gene"].unique()))
    candf = candf.merge(_A[["gene", "protein_name", "categories", "uniprot_function"]], on="gene", how="left")
    candf["categories"] = candf["categories"].fillna("")
    _out = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/raw_manhattan_newpeak_genes_{CLS}_site_tile_FDR.csv"
    candf.to_csv(_out, index=False)
    _ntag = int((candf["categories"].str.len() > 0).sum())
    print(f"\n{len(candf)} unique {CLS}-only FDR(q<{QCUT}) new-peak genes across 20 axes "
          f"({_ntag} climate/stress/flowering/defense-tagged); wrote {os.path.basename(_out)}")
    print("category counts:")
    print(candf["categories"].str.split(",").explode().replace("", pd.NA).dropna().value_counts().to_string())
    display(candf[["gene", "protein_name", "categories", "block", "chrom", "pos", "n_axes", "q", "nlp"]].head(80))
else:
    print("no gene-bearing FDR new peaks on any axis")
'''


def fdr_code(cls):
    return _FDR_CODE.replace("__CLS__", cls)


def append_fdr(cls, title_word):
    """Append the FDR section to an already-built+executed site notebook, executing ONLY the
    new cells (so the manhattans are NOT re-run)."""
    out = os.path.join(HERE, f"raw_manhattan_snp_vs_{cls}_lfmm_site_tile.ipynb")
    nb = nbf.read(out, as_version=4)
    tmp = new_notebook()
    tmp["cells"] = [new_markdown_cell(fdr_md(cls, title_word)), new_code_cell(fdr_code(cls))]
    ExecutePreprocessor(timeout=3600, kernel_name="basic", startup_timeout=180).preprocess(
        tmp, {"metadata": {"path": HERE}})
    nb["cells"].extend(tmp["cells"])
    with open(out, "w") as f:
        nbf.write(nb, f)
    print(f"[appended FDR section, manhattans untouched] {out}")


def build(cls, title_word):
    C = []
    md = lambda s: C.append(new_markdown_cell(s))
    co = lambda s: C.append(new_code_cell(s))

    md(f"""# Site-collapsed raw (before-WZA) mirror Manhattans — SNP vs {title_word}, LFMM, clq0.9 TILING, all 20 axes

**Site twin** of `raw_manhattan_snp_vs_{cls}_lfmm_tile.ipynb`. The per-record LFMM p here is
computed on the **site-mean Δp** — the mean allele-frequency change across the ~11 pools
within each of the **31 sites** — then `lfmm_ridge(K=1)` + `lfmm_test`. The pool notebook
treats all 352 pools as independent; that pseudoreplication inflates the pool scan
(genomic-inflation λ≈2.5). Collapsing pools→sites removes that design effect (site λ≈1.7–1.8)
and, more usefully, acts as a **filter**: a pool peak driven by a few correlated pools within
a site *collapses* at site level, while a genuinely site-coherent climate association *persists*.
`K=1` is the site-level operating point (at n=31, larger K re-inflates by overfitting).

Per-record `wza_in_clq09_tile_site` p-values (site-level RAW lfmm p, `pval`) — one figure per
climate axis (bio1-19 + pc1), SNP up, {cls} down. Dashed = per-class Bonferroni (0.05/n records).
QQ panel shows the per-axis genomic-inflation λ (compare to the pool notebook's ≈2–2.5).
Vertical dotted line + gene label = a tiling block whose lead {cls} record is Bonferroni-sig
while its lead SNP record in the same block is not ("new" peak, first-pass).""")

    co(f"""import os, sys
import numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
import lib
GEA = lib.GEA
WZAIN = f"{{GEA}}/r2_gea_nonsnp/phase1_replication/results/multiaxis/wza_in_clq09_tile_site"
AXES = [f"bio{{i}}" for i in range(1, 20)] + ["pc1"]
MODEL = "{MODEL}"
CLS = "{cls}"
CHROM_LEN = {{"Chr1":30427671,"Chr2":19698289,"Chr3":23459830,"Chr4":18585056,"Chr5":26975502}}
CH = list(CHROM_LEN); OFF = {{}}; _a = 0
for _c in CH: OFF[_c] = _a; _a += CHROM_LEN[_c]
plt.rcParams.update({{"figure.dpi":110}})
SNP_DARK, SNP_LIGHT = "#495057", "#ADB5BD"
CLS_DARK  = "#c0392b" if CLS == "sv" else "#1B5E20"
CLS_LIGHT = "#e8a58c" if CLS == "sv" else "#66BB6A"

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

    co("""from scipy import stats

def load_raw(cls, axis):
    f = f"{WZAIN}/{MODEL}_{cls}_gen9_{axis}.csv"
    if not os.path.exists(f): return None
    d = pd.read_csv(f)
    d = d[d["MAF"] > 0.05].copy()
    d["nlp"] = -np.log10(d["pval"].clip(lower=1e-300))
    d["gpos"] = d["chrom"].map(OFF) + d["pos"]
    return d

def leads(d):
    d = d[d["block"].notna() & (d["block"] != "")]
    return d.loc[d.groupby("block")["nlp"].idxmax()]

def gif_lambda(p):
    # genomic inflation factor: median observed chi2(1df) / expected median (0.4549)
    p = np.clip(np.asarray(p, float), 1e-300, 1.0)
    chi2 = stats.chi2.isf(p, 1)
    return float(np.median(chi2) / stats.chi2.isf(0.5, 1))

def qq_xy(p, keep_head=15000, n_bulk=3000):
    # observed vs expected -log10p under the uniform null; thin the dense bulk but
    # keep the whole significant tail (small p). Returns (expected, observed).
    p = np.sort(np.clip(np.asarray(p, float)[np.isfinite(p)], 1e-300, 1.0))
    n = len(p)
    exp = -np.log10((np.arange(n) + 0.5) / n)
    obs = -np.log10(p)
    if n <= keep_head + n_bulk:
        idx = np.arange(n)
    else:
        head = np.arange(keep_head)                                   # strongest hits, kept in full
        bulk = np.unique(np.geomspace(keep_head, n - 1, n_bulk).astype(int))
        idx = np.concatenate([head, bulk])
    return exp[idx], obs[idx]

def raw_mirror(ax, axis, cls, s, c):
    bonf_s, bonf_c = -np.log10(0.05/len(s)), -np.log10(0.05/len(c))
    ls_, lc_ = leads(s), leads(c)
    m = lc_.merge(ls_[["block", "nlp"]], on="block", suffixes=("", "_snp"), how="left")
    m["nlp_snp"] = m["nlp_snp"].fillna(0.0)
    newpk = m[(m["nlp"] > bonf_c) & (m["nlp_snp"] <= bonf_s)].copy()
    for i, ch in enumerate(CH):
        ss, cc = s[s.chrom == ch], c[c.chrom == ch]
        ax.scatter(ss.gpos, ss.nlp, s=2, c=(SNP_DARK if i % 2 == 0 else SNP_LIGHT), rasterized=True, linewidths=0)
        ax.scatter(cc.gpos, -cc.nlp, s=2, c=(CLS_DARK if i % 2 == 0 else CLS_LIGHT), rasterized=True, linewidths=0)
    ax.axhline(0, c="k", lw=.7)
    ax.axhline(bonf_s, ls="--", c="0.35", lw=.9)
    ax.axhline(-bonf_c, ls="--", c=CLS_DARK, lw=.9)
    ymax = max(s.nlp.max(), c.nlp.max()) * 1.15
    newpk["genes"] = newpk["block"].map(lambda b: ";".join(genes_on(b)))
    for _, r in newpk.iterrows():
        lab = r.genes.split(";")[0] if r.genes else r.block
        ax.axvline(r.gpos, ls=":", c="0.55", lw=.7)
        ax.annotate(lab, xy=(r.gpos, -ymax*0.95), ha="center", va="bottom", fontsize=6,
                    color=CLS_DARK, rotation=90)
    ax.set_ylim(-ymax, ymax)
    ax.set_xticks([OFF[c]+CHROM_LEN[c]/2 for c in CH]); ax.set_xticklabels(CH, fontsize=8)
    ax.set_ylabel(f"-log10 p  (SNP \\u2191   {cls} \\u2193)", fontsize=9)
    ax.annotate(f"{axis} \\u00b7 lfmm-SITE \\u00b7 clq0.9 tiling \\u00b7 {len(newpk)} {cls}-only Bonferroni new peaks",
                xy=(0.004, 0.97), xycoords="axes fraction", va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.7", alpha=.85))
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    return newpk[["block", "chrom", "pos", "nlp", "nlp_snp", "genes"]].sort_values("nlp", ascending=False)

def qq_panel(ax, cls, s, c):
    lam_s, lam_c = gif_lambda(s["pval"]), gif_lambda(c["pval"])
    xs, ys = qq_xy(s["pval"].to_numpy()); xc, yc = qq_xy(c["pval"].to_numpy())
    hi = max(xs.max(), ys.max(), xc.max(), yc.max())
    ax.plot([0, hi], [0, hi], c="k", lw=.8, ls="--", zorder=1)                 # null diagonal
    ax.scatter(xs, ys, s=4, c=SNP_DARK, rasterized=True, linewidths=0, zorder=2,
               label=f"SNP  (\\u03bb={lam_s:.2f})")
    ax.scatter(xc, yc, s=4, c=CLS_DARK, rasterized=True, linewidths=0, zorder=2,
               label=f"{cls}  (\\u03bb={lam_c:.2f})")
    ax.set_xlabel("expected -log10 p", fontsize=9)
    ax.set_ylabel("observed -log10 p", fontsize=9)
    ax.legend(loc="upper left", frameon=False, fontsize=8, markerscale=2, handletextpad=0.3)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    return lam_s, lam_c

all_newpk = []
for axis in AXES:
    s, c = load_raw("snp", axis), load_raw(CLS, axis)
    if s is None or c is None:
        print(f"missing {axis}"); continue
    fig = plt.figure(figsize=(19, 5.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.5, 1], wspace=0.18)
    newpk = raw_mirror(fig.add_subplot(gs[0]), axis, CLS, s, c)
    qq_panel(fig.add_subplot(gs[1]), CLS, s, c)
    fig.tight_layout(); plt.show()
    newpk = newpk.copy(); newpk["axis"] = axis
    all_newpk.append(newpk)
    print(f"{axis}: {len(newpk)} {CLS}-only Bonferroni new peaks (gene-bearing: {(newpk.genes != '').sum()})")
    if len(newpk):
        display(newpk[["block", "chrom", "pos", "genes", "nlp", "nlp_snp"]].rename(
            columns={"nlp": f"{CLS}_nlp", "nlp_snp": "snp_nlp_here"}))""")

    md(f"""## Candidate genes — {title_word}-only new peaks (SITE-collapsed), consolidated across all 20 axes

Union of every gene-bearing new peak found above, deduplicated per gene (kept: the
axis with the strongest {cls} signal), with Ensembl Plants symbol/description and a
recurrence count (how many of the 20 axes flag that gene's block). **Site-collapsed raw p**
(pseudoreplication removed by pooling to 31 sites; no GIF applied to the plotted p). Compare
this list to the pool `raw_manhattan_newpeak_genes_{cls}_tile.csv` — genes that drop out here
were pool-pseudoreplication artifacts; genes that persist are site-coherent candidates.""")

    co('''import json, time, urllib.request

ALL = pd.concat(all_newpk, ignore_index=True) if all_newpk else pd.DataFrame(columns=["block","gene","axis","nlp"])
gene_rows = []
for _, r in ALL.iterrows():
    for g in (r.genes.split(";") if r.genes else []):
        gene_rows.append(dict(gene=g, block=r.block, chrom=r.chrom, pos=r.pos, axis=r.axis, nlp=r.nlp))
G = pd.DataFrame(gene_rows)
if len(G):
    recur = G.groupby("gene")["axis"].nunique().rename("n_axes")
    best = G.sort_values("nlp", ascending=False).drop_duplicates("gene").set_index("gene")
    cand = best.join(recur).sort_values("n_axes", ascending=False).reset_index()

    def ensembl_symbols(ids, chunk=900, retries=3):
        out = {}
        for i in range(0, len(ids), chunk):
            sub = ids[i:i + chunk]
            req = urllib.request.Request("https://rest.ensembl.org/lookup/id",
                data=json.dumps({"ids": sub}).encode(),
                headers={"Content-Type": "application/json", "Accept": "application/json"})
            for attempt in range(retries):
                try:
                    with urllib.request.urlopen(req, timeout=60) as r:
                        d = json.load(r)
                    for k, v in d.items():
                        if v:
                            out[k] = (v.get("display_name", "") or "", (v.get("description", "") or "").split(" [Source")[0])
                    break
                except Exception as e:
                    print(f"  [ensembl] chunk {i} attempt {attempt+1}: {e}", flush=True); time.sleep(2)
        return out

    sym = ensembl_symbols(sorted(cand.gene.unique()))
    cand["symbol"] = cand.gene.map(lambda x: sym.get(x, ("", ""))[0])
    cand["description"] = cand.gene.map(lambda x: sym.get(x, ("", ""))[1])
    cand.to_csv(f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/raw_manhattan_newpeak_genes_{CLS}_site_tile.csv", index=False)
    print(f"{len(cand)} unique {CLS}-only new-peak genes (site-collapsed) across 20 axes")
    display(cand[["gene", "symbol", "description", "block", "chrom", "pos", "n_axes", "nlp"]].head(50))
else:
    print("no gene-bearing new peaks found on any axis")''')

    md(f"""## Functional theme per gene — stress / climate / flowering / circadian / defense

Same {title_word}-only new-peak gene list (site-collapsed), tagged with a **functional
category** via the project annotator `annotate_genes_tair_uniprot.py` (TAIR GO + UniProt —
the mandated route). Categories: `flowering` (flowering time, vernalization, photoperiod,
**circadian**, meristem identity) · `temperature` (cold/freezing/heat, CBF/DREB/COR) ·
`water` (drought, osmotic, ABA, salinity, stomata) · `light` (photomorphogenesis, UV, shade
avoidance, photosynthesis) · `oxidative` (ROS, abiotic/general stress) · `defense` (immunity,
pathogen, SA/JA, R-genes) · `calcium` (Ca2+/calmodulin/CDPK). A gene can carry several.
Written to `raw_manhattan_newpeak_genes_{cls}_site_tile_annotated.csv`.""")

    co('''import importlib.util as _ilu

_annp = os.path.abspath(os.path.join(os.getcwd(), "..", "annotate_genes_tair_uniprot.py"))
_sp = _ilu.spec_from_file_location("ann_tu", _annp)
ann_tu = _ilu.module_from_spec(_sp); _sp.loader.exec_module(ann_tu)

if "cand" in dir() and len(cand):
    _A = ann_tu.annotate(sorted(cand["gene"].unique()))
    _cols = ["gene", "protein_name", "categories", "climate_stress_flowering",
             "uniprot_function", "uniprot_keywords"]
    cand_annot = cand.merge(_A[_cols], on="gene", how="left")
    cand_annot["categories"] = cand_annot["categories"].fillna("")
    _outp = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/raw_manhattan_newpeak_genes_{CLS}_site_tile_annotated.csv"
    cand_annot.to_csv(_outp, index=False)
    _ntag = int((cand_annot["categories"].str.len() > 0).sum())
    print(f"{len(cand_annot)} genes annotated ; {_ntag} tagged to >=1 climate/stress/flowering/defense category")
    print("\\ncategory counts (a gene can carry several):")
    print(cand_annot["categories"].str.split(",").explode().replace("", pd.NA)
          .dropna().value_counts().to_string())
    display(cand_annot[["gene", "symbol", "protein_name", "categories", "n_axes",
                        "nlp", "chrom", "pos"]].sort_values(["n_axes", "nlp"],
                        ascending=False).head(60))
else:
    cand_annot = pd.DataFrame(columns=["gene", "symbol", "categories"])
    print("no genes to annotate")''')

    md(f"""## Dot plot — recurrent {title_word}-only new-peak genes (site-collapsed) across the 20 climate axes

One dot per (gene, axis) where that gene's block is an {cls}-only Bonferroni new peak on
that axis. **Dot size** = -log10 p (that axis' {cls} lead); **color** = functional category
(above). Genes are ordered genomically (chrom then position) and the panel is limited to
genes recurring on **enough axes to stay legible** (auto-threshold, printed below); the full
gene list with categories is the table above. Reading across a row shows how many of the 20
climate variables independently flag that gene's block at site level.""")

    co('''import matplotlib.lines as mlines

AXES_ORD = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
_rows = []
if all_newpk:
    for r in pd.concat(all_newpk, ignore_index=True).itertuples():
        for g in (r.genes.split(";") if r.genes else []):
            _rows.append(dict(gene=g, axis=r.axis, nlp=r.nlp, chrom=r.chrom, pos=r.pos))
GX = pd.DataFrame(_rows)

if len(GX):
    _nax = GX.groupby("gene")["axis"].nunique()
    thr = 2
    while int((_nax >= thr).sum()) > 70 and thr < 8:
        thr += 1
    keep_genes = set(_nax[_nax >= thr].index)
    D = GX[GX["gene"].isin(keep_genes)].drop_duplicates(["gene", "axis"]).copy()
    gord = (D.drop_duplicates("gene")
            .assign(_cn=lambda x: x["chrom"].str.replace("Chr", "", regex=False).astype(int))
            .sort_values(["_cn", "pos"])["gene"].tolist())
    catmap = dict(zip(cand_annot["gene"], cand_annot["categories"])) if len(cand_annot) else {}
    symmap = dict(zip(cand_annot["gene"], cand_annot["symbol"])) if len(cand_annot) else {}
    def _pcat(g):
        c = [x for x in (catmap.get(g, "") or "").split(",") if x]
        return c[0] if c else "unclassified"
    CATCOL = {"flowering": "#8E44AD", "temperature": "#2980B9", "water": "#16A085",
              "light": "#F1C40F", "oxidative": "#E67E22", "defense": "#C0392B",
              "calcium": "#7F8C8D", "unclassified": "#CED4DA"}
    yidx = {g: i for i, g in enumerate(gord)}
    xidx = {a: i for i, a in enumerate(AXES_ORD)}
    fig, ax = plt.subplots(figsize=(9, max(4, 0.26 * len(gord))))
    for _, rr in D.iterrows():
        ax.scatter(xidx[rr["axis"]], yidx[rr["gene"]], s=12 + rr["nlp"] * 2.2,
                   c=CATCOL.get(_pcat(rr["gene"]), "#CED4DA"), edgecolors="white",
                   linewidths=0.3, alpha=0.95, zorder=3)
    ax.set_xticks(range(len(AXES_ORD))); ax.set_xticklabels(AXES_ORD, rotation=90, fontsize=7)
    ax.set_yticks(range(len(gord)))
    ax.set_yticklabels([f"{(symmap.get(g) or g)}  ({_pcat(g)})" for g in gord], fontsize=6)
    ax.set_xlim(-0.6, len(AXES_ORD) - 0.4); ax.set_ylim(-0.7, len(gord) - 0.3)
    ax.grid(True, axis="both", lw=0.3, c="0.9", zorder=0); ax.set_axisbelow(True)
    ax.set_xlabel("climate axis", fontsize=9)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    _present = [c for c in CATCOL if c in {_pcat(g) for g in gord}]
    _h = [mlines.Line2D([], [], marker="o", ls="", ms=7, mfc=CATCOL[c], mec="white", label=c)
          for c in _present]
    ax.legend(handles=_h, title="category", loc="upper left", bbox_to_anchor=(1.01, 1.0),
              frameon=False, fontsize=7, title_fontsize=8)
    fig.tight_layout(); plt.show()
    _tot = GX["gene"].nunique()
    print(f"dot plot: {len(gord)} genes recurring on >={thr} of 20 axes "
          f"(auto-threshold for legibility); {_tot - len(gord)} of {_tot} lower-recurrence "
          f"genes omitted from the plot but present in the annotated table/CSV.")
else:
    print("no gene-bearing new peaks to plot")''')

    md(fdr_md(cls, title_word)); co(fdr_code(cls))     # FDR-relaxed gene list at the end

    nb = new_notebook(); nb["cells"] = C
    out = os.path.join(HERE, f"raw_manhattan_snp_vs_{cls}_lfmm_site_tile.ipynb")
    ep = ExecutePreprocessor(timeout=3600, kernel_name="basic", startup_timeout=180)
    ep.preprocess(nb, {"metadata": {"path": HERE}})
    with open(out, "w") as f:
        nbf.write(nb, f)
    print(f"[built+executed] {out}")


if __name__ == "__main__":
    import sys
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("all", "nonsnp"):
        build("nonsnp", "non-SNP (pooled)")   # full build (manhattans + Bonferroni genes + FDR genes)
    if mode in ("all", "sv_fdr"):
        append_fdr("sv", "SV")                 # append FDR section to existing SV nb (manhattans untouched)
