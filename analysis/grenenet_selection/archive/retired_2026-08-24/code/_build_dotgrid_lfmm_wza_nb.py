#!/usr/bin/env python
"""WZA (block-level) twin of the class-peak x climate-axis dot grid.

Rows = class blocks that are WZA-Bonferroni-sig on >=1 axis (isotonic WZA on the no-GIF
LFMM p). Columns = 20 climate axes. Per (block, axis):
  - dot size = block -log10(Z_pVal) for the class
  - BLACK RING = block WZA-Bonf-sig in the class but NOT in the co-located SNP block
  - colour = per-block max LD r2 of that block's class anchors to a SNP within +/-50 kb
    (masked `sv_snp_ld_v2_maconly`; uniform down a row since it is a block property;
    dark = the block's variants are tightly tagged by a nearby SNP; pale = not)
Two figures: all peaks (dense, capped for readability) + top-45 compact. Emits both:
  sv_peak_dotgrid_lfmm_wza.ipynb   (rows = SV blocks; r2 over SV anchors)
  nonsnp_peak_dotgrid_lfmm_wza.ipynb (rows = non-SNP blocks; r2 over SV+indel anchors)
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def build(cls, title_word, r2_kinds, dense_cap):
    C = []
    md = lambda s: C.append(new_markdown_cell(s))
    co = lambda s: C.append(new_code_cell(s))

    md(f"""# {title_word} peaks x climate axes — WZA dot grid (SNP-tagging view)

Block-level twin of `{cls}_peak_dotgrid_lfmm.ipynb`, using the **WZA** (isotonic) block
result on the no-GIF LFMM p. **Rows** = {title_word} blocks WZA-Bonferroni-sig on >=1
axis. **Columns** = 20 climate axes. **Dot size** = block -log10(Z_pVal); **black ring** =
block WZA-Bonf-sig in {title_word} but NOT in the co-located SNP block; **colour = per-block
max LD r2** of the block's {title_word} anchors to a SNP within +/-50 kb (masked
`sv_snp_ld_v2_maconly`; uniform down each row — it is a block property). Pale ring =
{title_word}-only WZA block whose variants are not tightly SNP-tagged.""")

    co(f'''import os, sys
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), "..", "..")))
import lib
GEA = lib.GEA
WD = f"{{GEA}}/phase1_replication/results/multiaxis/wza"
OUTDIR = f"{{GEA}}/phase1_replication/results/multiaxis"
AXES = [f"bio{{i}}" for i in range(1, 20)] + ["pc1"]
CLS = "{cls}"; R2_KINDS = {r2_kinds!r}; DENSE_CAP = {dense_cap}
plt.rcParams.update({{"figure.dpi": 200}})

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
print("loaded:", len(SPANS), "blocks,", len(GENES), "genes ; class =", CLS)''')

    co('''# --- per-anchor tagging r2 (masked, no-missingness variant), for the block max-r2 ---
LDDIR = f"{GEA}/sv_snp_ld_v2_maconly"
anc = {ci: {"pos": [], "r2": []} for ci in range(1, 6)}
for kind in R2_KINDS:
    for ci in range(1, 6):
        f = f"{LDDIR}/tagging_{kind}_panel_Chr{ci}.npz"
        if not os.path.exists(f): continue
        z = np.load(f, allow_pickle=True)
        anc[ci]["pos"].append(z["pos"].astype(np.int64)); anc[ci]["r2"].append(z["best_r2"].astype(float))
ANC = {}
for ci in range(1, 6):
    if anc[ci]["pos"]:
        p = np.concatenate(anc[ci]["pos"]); r = np.concatenate(anc[ci]["r2"])
        o = np.argsort(p); ANC[f"Chr{ci}"] = (p[o], r[o])
def block_r2(block):
    if block not in SPANS.index: return np.nan
    s = SPANS.loc[block]; key = s.chrom
    if key not in ANC: return np.nan
    P, R = ANC[key]; lo, hi = np.searchsorted(P, s.start), np.searchsorted(P, s.end, side="right")
    seg = R[lo:hi]; seg = seg[np.isfinite(seg)]
    return float(seg.max()) if len(seg) else np.nan

def load_wza(cls, axis):
    f = f"{WD}/wza_lfmm_{cls}_gen9_{axis}_isotonic.csv"
    if not os.path.exists(f): return None
    w = pd.read_csv(f).rename(columns={"index": "block"})
    w["block"] = w["block"].astype(str)
    w = w[w["Z_pVal"].notna() & (w["Z_pVal"] > 0)].copy()
    w["nlp"] = -np.log10(w["Z_pVal"].clip(lower=1e-300))
    return w

# grid[(block, axis)] = (nlp, cls_sig, snp_sig)
grid = {}; recur = {}
for axis in AXES:
    c = load_wza(CLS, axis); s = load_wza("snp", axis)
    if c is None or s is None: continue
    bonf_c = -np.log10(0.05/len(c)); bonf_s = -np.log10(0.05/len(s))
    ssig = set(s.loc[s["nlp"] > bonf_s, "block"])
    for _, row in c.iterrows():
        b = row["block"]; v = float(row["nlp"]); cls_sig = v > bonf_c; snp_sig = b in ssig
        grid[(b, axis)] = (v, cls_sig, snp_sig)
        if cls_sig: recur[b] = recur.get(b, 0) + 1
peaks = sorted(recur, key=lambda b: (recur[b], b), reverse=True)
R2B = {b: block_r2(b) for b in peaks}
print(f"{len(peaks)} {CLS} WZA-Bonferroni blocks on >=1 of 20 axes")''')

    co('''def rlabel(b):
    gs = genes_on(b); s = SPANS.loc[b] if b in SPANS.index else None
    pos = f"{s.chrom}:{s.start/1e6:.1f}" if s is not None else b
    return f"{gs[0] if gs else b} {pos}"

norm = mpl.colors.Normalize(vmin=0, vmax=1.0); cmap = plt.cm.OrRd
NA_COLOR = "0.82"

def draw(peaks_sub, pitch, dot_base, dot_scale, yfont, ringlw):
    if not peaks_sub:
        print("(no blocks to plot)"); return
    fig, ax = plt.subplots(figsize=(9.0, 0.7 + pitch * len(peaks_sub)))
    n_ring = n_na = 0
    for xi, axis in enumerate(AXES):
        for yi, b in enumerate(peaks_sub):
            g = grid.get((b, axis))
            if g is None: continue
            nlp, cls_sig, snp_sig = g
            ring = cls_sig and not snp_sig; n_ring += int(ring)
            r2 = R2B.get(b, np.nan)
            col = NA_COLOR if not np.isfinite(r2) else cmap(norm(r2)); n_na += int(not np.isfinite(r2))
            ax.scatter(xi, yi, s=dot_base + dot_scale * min(nlp, 10), c=[col],
                       edgecolors=("black" if ring else "none"), linewidths=(ringlw if ring else 0), zorder=3)
    ax.set_xticks(range(len(AXES))); ax.set_xticklabels(AXES, rotation=90, fontsize=5)
    ax.set_yticks(range(len(peaks_sub))); ax.set_yticklabels([rlabel(b) for b in peaks_sub], fontsize=yfont)
    ax.tick_params(length=1.5, pad=1.0)
    ax.set_xlabel("climate axis", fontsize=7); ax.set_ylim(-0.6, len(peaks_sub) - 0.4); ax.margins(x=0.02)
    ax.grid(True, alpha=0.13, lw=0.3, zorder=0)
    ax.annotate(f"lfmm+WZA \\u00b7 {CLS} \\u2014 size \\u221d block \\u2212log10p    \\u25cb = block not shared by a SNP peak "
                f"({n_ring} ringed / {len(peaks_sub)} blocks)", (0.0, 1.004), xycoords="axes fraction", fontsize=5.5, va="bottom")
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.03, pad=0.01)
    cb.set_label("per-block max LD r\\u00b2 to a SNP within 50 kb\\n(dark = block variants tagged by a nearby SNP)", fontsize=6)
    cb.ax.tick_params(labelsize=5, length=2)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    fig.tight_layout(); plt.show()
    print(f"{len(peaks_sub)} blocks: {n_ring} ringed cells; {n_na} grey (untestable)")''')

    md(f"""## All {title_word} WZA blocks (dense)""")
    co("""if len(peaks) > DENSE_CAP:
    print(f'dense view: showing {DENSE_CAP} of {len(peaks)} blocks (full set in CSV)')
    draw(peaks[:DENSE_CAP], pitch=0.125, dot_base=2, dot_scale=7, yfont=3.0, ringlw=0.55)
else:
    draw(peaks, pitch=0.16, dot_base=3, dot_scale=9, yfont=4.0, ringlw=0.7)""")

    md(f"""## Top 45 most-recurrent {title_word} WZA blocks (compact)""")
    co("draw(peaks[:45], pitch=0.20, dot_base=3, dot_scale=11, yfont=4.5, ringlw=0.8)")

    co(f'''rows = []
for (b, axis), (nlp, cls_sig, snp_sig) in grid.items():
    if b not in peaks: continue
    gs = genes_on(b)
    rows.append(dict(block=b, axis=axis, gene=(gs[0] if gs else ""), block_nlp=round(nlp,3),
                     cls_bonf_sig=cls_sig, snp_bonf_sig=snp_sig, cls_only=(cls_sig and not snp_sig),
                     block_max_snp_r2=(round(R2B.get(b,np.nan),4) if np.isfinite(R2B.get(b,np.nan)) else np.nan)))
tab = pd.DataFrame(rows).sort_values(["block","axis"]) if rows else pd.DataFrame()
tab.to_csv(f"{{OUTDIR}}/{cls}_peak_dotgrid_lfmm_wza.csv", index=False)
print("wrote", f"{{OUTDIR}}/{cls}_peak_dotgrid_lfmm_wza.csv", "|", len(peaks), "blocks x", len(AXES), "axes")''')

    nb = new_notebook(); nb["cells"] = C
    out = os.path.join(HERE, f"{cls}_peak_dotgrid_lfmm_wza.ipynb")
    ep = ExecutePreprocessor(timeout=3600, kernel_name="basic", startup_timeout=180)
    ep.preprocess(nb, {"metadata": {"path": HERE}})
    with open(out, "w") as f:
        nbf.write(nb, f)
    print(f"[built+executed] {out}")


if __name__ == "__main__":
    build("sv", "SV", ["sv"], dense_cap=150)
    build("nonsnp", "non-SNP", ["sv", "indel"], dense_cap=150)
