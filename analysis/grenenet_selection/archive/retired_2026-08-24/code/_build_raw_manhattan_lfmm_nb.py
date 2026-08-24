#!/usr/bin/env python
"""Build raw (BEFORE-WZA) mirror Manhattans, LFMM only, all 20 climate axes.

Per-record wza_in pvals (not the WZA block aggregate) -- the same "SNP up / class
down" style as the per-garden persite notebook (dense per-variant scatter, dashed
Bonferroni lines, vertical dotted lines + gene labels at class-only Bonferroni new
peaks). "New peak" here = a clq0.9 block whose LEAD (min-p) record in the class is
Bonferroni-sig while its lead SNP record in the same block is not.

Emits two twin notebooks (same code, class swapped):
  raw_manhattan_snp_vs_nonsnp_lfmm.ipynb  -- SNP vs pooled non-SNP
  raw_manhattan_snp_vs_sv_lfmm.ipynb      -- SNP vs SV
"""
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

    md(f"""# Raw (before-WZA) mirror Manhattans — SNP vs {title_word}, LFMM, all 20 axes

Per-record `wza_in` p-values (the input WZA aggregates, not the block-level output) --
one figure per climate axis (bio1-19 + pc1), LFMM only. SNP up, {cls} down. Dashed =
per-class Bonferroni (0.05/n records). Vertical dotted line + gene label = a clq0.9
block whose lead {cls} record is Bonferroni-sig while its lead SNP record in the same
block is not ("new" peak, first-pass -- not yet checked for SNP-tagging).""")

    co(f"""import os, sys
import numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), "..", "..")))
import lib
GEA = lib.GEA
WZAIN = f"{{GEA}}/phase1_replication/results/multiaxis/wza_in"
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

def block_spans(r2=0.9):
    tag=f"clq{{r2}}"; rows=[]
    for ci in range(1,6):
        g=pd.read_csv(f"{{lib.CLQ_BLOCKS_DIR}}/chr{{ci}}_{{tag}}_blocks_{{tag}}.tsv",sep="\\t").sort_values("start_pos").reset_index(drop=True)
        for idx,r in g.iterrows():
            rows.append((f"Chr{{ci}}_{{idx}}",f"Chr{{ci}}",int(r.start_pos),int(r.end_pos)))
    return pd.DataFrame(rows,columns=["block","chrom","start","end"]).set_index("block")
SPANS=block_spans(); GENES=lib.load_genes()
def genes_on(block):
    if block not in SPANS.index: return []
    s=SPANS.loc[block]; gc=GENES[(GENES.chrom==s.chrom)&(GENES.end>=s.start)&(GENES.start<=s.end)]
    return list(gc.gene)
print("loaded:", len(SPANS), "blocks,", len(GENES), "genes ; class =", CLS)""")

    co("""from scipy import stats

def load_raw(cls, axis):
    f = f"{WZAIN}/{MODEL}_{cls}_gen9_{axis}.csv"
    if not os.path.exists(f): return None
    d = pd.read_csv(f)
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
    newpk = m[(m["nlp"] > bonf_c) & (m["nlp_snp"] <= bonf_s)]
    for i, ch in enumerate(CH):
        ss, cc = s[s.chrom == ch], c[c.chrom == ch]
        ax.scatter(ss.gpos, ss.nlp, s=2, c=(SNP_DARK if i % 2 == 0 else SNP_LIGHT), rasterized=True, linewidths=0)
        ax.scatter(cc.gpos, -cc.nlp, s=2, c=(CLS_DARK if i % 2 == 0 else CLS_LIGHT), rasterized=True, linewidths=0)
    ax.axhline(0, c="k", lw=.7)
    ax.axhline(bonf_s, ls="--", c="0.35", lw=.9)
    ax.axhline(-bonf_c, ls="--", c=CLS_DARK, lw=.9)
    ymax = max(s.nlp.max(), c.nlp.max()) * 1.15
    for _, r in newpk.iterrows():
        gs = genes_on(r.block); lab = gs[0] if gs else r.block
        ax.axvline(r.gpos, ls=":", c="0.55", lw=.7)
        ax.annotate(lab, xy=(r.gpos, -ymax*0.95), ha="center", va="bottom", fontsize=6,
                    color=CLS_DARK, rotation=90)
    ax.set_ylim(-ymax, ymax)
    ax.set_xticks([OFF[c]+CHROM_LEN[c]/2 for c in CH]); ax.set_xticklabels(CH, fontsize=8)
    ax.set_ylabel(f"-log10 p  (SNP \\u2191   {cls} \\u2193)", fontsize=9)
    ax.annotate(f"{axis} \\u00b7 lfmm \\u00b7 {len(newpk)} {cls}-only Bonferroni new peaks",
                xy=(0.004, 0.97), xycoords="axes fraction", va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.7", alpha=.85))
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    return len(newpk)

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

for axis in AXES:
    s, c = load_raw("snp", axis), load_raw(CLS, axis)
    if s is None or c is None:
        print(f"missing {axis}"); continue
    fig = plt.figure(figsize=(19, 5.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.5, 1], wspace=0.18)
    raw_mirror(fig.add_subplot(gs[0]), axis, CLS, s, c)
    qq_panel(fig.add_subplot(gs[1]), CLS, s, c)
    fig.tight_layout(); plt.show()""")

    nb = new_notebook(); nb["cells"] = C
    out = os.path.join(HERE, f"raw_manhattan_snp_vs_{cls}_lfmm.ipynb")
    ep = ExecutePreprocessor(timeout=3600, kernel_name="basic", startup_timeout=180)
    ep.preprocess(nb, {"metadata": {"path": HERE}})
    with open(out, "w") as f:
        nbf.write(nb, f)
    print(f"[built+executed] {out}")


if __name__ == "__main__":
    build("nonsnp", "non-SNP (pooled)")
    build("sv", "SV")
