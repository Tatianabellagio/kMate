#!/usr/bin/env python
"""GIF-corrected twin of the non-SNP-peak x climate-axis dot grid (nonsnp_peak_dotgrid_lfmm.ipynb).

Identical to the raw dot grid EXCEPT the per-record LFMM p-values are CONDITIONALLY
genomic-inflation-corrected before peaks / significance / dot sizes are computed:
    lambda = median chi2(1df) / 0.4549   (per axis, per class)
    lambda > 1  -> chi2/lambda, recompute p (deflate the inflated axes)
    lambda <= 1 -> leave raw (correcting would inflate, e.g. bio19 ~0.92)
So an SV "peak" here = lead non-SNP record Bonferroni-sig AFTER this correction. Colour still =
max LD r2 of the variant to a SNP within 50 kb (masked, no-missingness `sv_snp_ld_v2_maconly`),
ring still = non-SNP peak not shared by a SNP peak. Two figures: all peaks (dense) + top-45.
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
import os

HERE = os.path.dirname(os.path.abspath(__file__))
C = []
md = lambda s: C.append(new_markdown_cell(s))
co = lambda s: C.append(new_code_cell(s))

md("""# non-SNP peaks x climate axes — GIF-corrected dot grid (SNP-tagging view)

GIF-corrected twin of `nonsnp_peak_dotgrid_lfmm.ipynb`. Same raw LFMM K=16 p-values, but
**conditionally genomic-inflation-corrected**: per axis x class, deflate by lambda when
lambda>1, leave raw when lambda<=1 (correcting a lambda<1 axis like bio19 would inflate).
**Rows** = non-SNP peak blocks (lead non-SNP record Bonferroni-sig *after correction* on >=1 axis).
**Columns** = 20 climate axes. **Dot size** = corrected SV -log10p; **colour = max LD r2 of
the lead non-SNP to a SNP within +/-50 kb** (masked `sv_snp_ld_v2_maconly`: dark = tagged by a
nearby SNP; pale = not); **black ring** = non-SNP peak not shared by a SNP peak (SV Bonf-sig,
no SNP in the block Bonf-sig). Far fewer peaks survive here than in the raw grid — that is
the point: most raw SV peaks were genomic inflation.""")

co('''import os, sys
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy import stats
sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), "..", "..")))
import lib
GEA = lib.GEA
WZAIN = f"{GEA}/phase1_replication/results/multiaxis/wza_in"
OUTDIR = f"{GEA}/phase1_replication/results/multiaxis"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
plt.rcParams.update({"figure.dpi": 200})
_EXP_MED = stats.chi2.isf(0.5, 1)

def block_spans(r2=0.9):
    tag=f"clq{r2}"; rows=[]
    for ci in range(1,6):
        g=pd.read_csv(f"{lib.CLQ_BLOCKS_DIR}/chr{ci}_{tag}_blocks_{tag}.tsv",sep="\\t").sort_values("start_pos").reset_index(drop=True)
        for idx,r in g.iterrows():
            rows.append((f"Chr{ci}_{idx}",f"Chr{ci}",int(r.start_pos),int(r.end_pos)))
    return pd.DataFrame(rows,columns=["block","chrom","start","end"]).set_index("block")
SPANS=block_spans(); GENES=lib.load_genes()
def genes_on(block):
    if block not in SPANS.index: return []
    s=SPANS.loc[block]; gc=GENES[(GENES.chrom==s.chrom)&(GENES.end>=s.start)&(GENES.start<=s.end)]
    return list(gc.gene)
print("loaded:", len(SPANS), "blocks,", len(GENES), "genes")''')

co('''# --- SV -> nearby-SNP LD tagging (best r2 within +/-50 kb), masked no-missingness variant ---
LDDIR = f"{GEA}/sv_snp_ld_v2_maconly"   # MAC floor >=2 kept, missingness gate dropped (min-pair-frac 0.0)
R2 = {}; POS2 = {}
for kind in ("sv", "indel"):
    for ci in range(1, 6):
        f = f"{LDDIR}/tagging_{kind}_panel_Chr{ci}.npz"
        if not os.path.exists(f): continue
        z = np.load(f, allow_pickle=True)
        ch = np.array([str(c) for c in z["chrom"]]); pos = z["pos"].astype(np.int64)
        size = z["size"].astype(np.int64); r2 = z["best_r2"].astype(float)
        for c, p, s, r in zip(ch, pos, size, r2):
            R2[(c, int(p), int(s))] = r
            if np.isfinite(r):
                k = (c, int(p)); POS2[k] = max(POS2.get(k, -1.0), r)
def lookup_r2(chrom, pos, size):
    r = R2.get((chrom, int(pos), int(size)))
    if r is not None and np.isfinite(r): return float(r)
    return POS2.get((chrom, int(pos)), np.nan)
print(f"tagging r2: {len(R2)} anchor records, {len(POS2)} testable positions")

def gif_correct(pval):
    p = np.clip(np.asarray(pval, float), 1e-300, 1.0)
    chi2 = stats.chi2.isf(p, 1); lam = float(np.median(chi2) / _EXP_MED)
    if lam > 1.0:
        return np.clip(stats.chi2.sf(chi2 / lam, 1), 1e-300, 1.0), lam, True
    return p, lam, False

def snp_lead_by_block(axis):
    d = pd.read_csv(f"{WZAIN}/lfmm_snp_gen9_{axis}.csv", usecols=["block", "pval"])
    d = d[d["block"].notna() & (d["block"] != "")]
    pc, lam, applied = gif_correct(d["pval"].to_numpy())
    d = d.assign(nlp=-np.log10(pc)); bonf = -np.log10(0.05 / len(d))
    return d.groupby("block")["nlp"].max(), bonf

def sv_lead_by_block(axis):
    d = pd.read_csv(f"{WZAIN}/lfmm_nonsnp_gen9_{axis}.csv",
                    usecols=["chrom", "pos", "ref_len", "alt_len", "block", "pval"])
    d = d[d["block"].notna() & (d["block"] != "")]
    pc, lam, applied = gif_correct(d["pval"].to_numpy())
    d = d.assign(nlp=-np.log10(pc)); bonf = -np.log10(0.05 / len(d))
    lead = d.loc[d.groupby("block")["nlp"].idxmax()].copy()
    lead["size"] = (lead["alt_len"] - lead["ref_len"]).abs()
    return lead, bonf

# grid[(block, axis)] = (sv_nlp_corr, sv_sig, snp_sig, r2)
grid = {}; recur = {}
for axis in AXES:
    sv_lead, sv_bonf = sv_lead_by_block(axis)
    snp_lead, snp_bonf = snp_lead_by_block(axis)
    for _, row in sv_lead.iterrows():
        b = row["block"]; v = float(row["nlp"])
        sn = float(snp_lead.get(b, 0.0)); sv_sig = v > sv_bonf; snp_sig = sn > snp_bonf
        r2 = lookup_r2(row["chrom"], row["pos"], row["size"])
        grid[(b, axis)] = (v, sv_sig, snp_sig, r2)
        if sv_sig: recur[b] = recur.get(b, 0) + 1
peaks = sorted(recur, key=lambda b: (recur[b], b), reverse=True)
print(f"{len(peaks)} non-SNP peak blocks (lead non-SNP Bonferroni-sig AFTER GIF on >=1 of 20 axes)")''')

co('''def rlabel(b):
    gs = genes_on(b); s = SPANS.loc[b] if b in SPANS.index else None
    pos = f"{s.chrom}:{s.start/1e6:.1f}" if s is not None else b
    return f"{gs[0] if gs else b} {pos}"

norm = mpl.colors.Normalize(vmin=0, vmax=1.0); cmap = plt.cm.OrRd
NA_COLOR = "0.82"

def draw(peaks_sub, pitch, dot_base, dot_scale, yfont, ringlw):
    if not peaks_sub:
        print("(no SV peaks survive GIF correction to plot)"); return
    fig, ax = plt.subplots(figsize=(9.0, 0.7 + pitch * len(peaks_sub)))
    n_ring = 0
    for xi, axis in enumerate(AXES):
        for yi, b in enumerate(peaks_sub):
            g = grid.get((b, axis))
            if g is None: continue
            sv_nlp, sv_sig, snp_sig, r2 = g
            ring = sv_sig and not snp_sig; n_ring += int(ring)
            col = NA_COLOR if not np.isfinite(r2) else cmap(norm(r2))
            ax.scatter(xi, yi, s=dot_base + dot_scale * min(sv_nlp, 10), c=[col],
                       edgecolors=("black" if ring else "none"), linewidths=(ringlw if ring else 0), zorder=3)
    ax.set_xticks(range(len(AXES))); ax.set_xticklabels(AXES, rotation=90, fontsize=5)
    ax.set_yticks(range(len(peaks_sub))); ax.set_yticklabels([rlabel(b) for b in peaks_sub], fontsize=yfont)
    ax.tick_params(length=1.5, pad=1.0)
    ax.set_xlabel("climate axis", fontsize=7); ax.set_ylim(-0.6, len(peaks_sub) - 0.4); ax.margins(x=0.02)
    ax.grid(True, alpha=0.13, lw=0.3, zorder=0)
    ax.annotate(f"GIF-corrected LFMM \\u00b7 non-SNP \\u2014 size \\u221d non-SNP \\u2212log10p    \\u25cb = non-SNP peak not shared by a SNP peak "
                f"({n_ring} ringed / {len(peaks_sub)} peaks)", (0.0, 1.004), xycoords="axes fraction", fontsize=5.5, va="bottom")
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.03, pad=0.01)
    cb.set_label("max LD r\\u00b2 of the variant to a SNP within 50 kb\\n(dark = tagged by a nearby SNP)", fontsize=6)
    cb.ax.tick_params(labelsize=5, length=2)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    fig.tight_layout(); plt.show()
    print(f"{len(peaks_sub)} peaks: {n_ring} ringed cells")''')

md("""## All non-SNP peaks surviving GIF correction (dense)""")
co("draw(peaks, pitch=0.125, dot_base=2, dot_scale=7, yfont=3.0, ringlw=0.55)")

md("""## Top 45 most-recurrent (compact)

If fewer than 45 survive, all of them are shown.""")
co("draw(peaks[:45], pitch=0.20, dot_base=3, dot_scale=11, yfont=4.5, ringlw=0.8)")

co('''rows = []
for (b, axis), (sv_nlp, sv_sig, snp_sig, r2) in grid.items():
    if b not in peaks: continue
    gs = genes_on(b)
    rows.append(dict(block=b, axis=axis, gene=(gs[0] if gs else ""), sv_nlp_gifcorr=round(sv_nlp,3),
                     sv_bonf_sig=sv_sig, snp_bonf_sig=snp_sig, sv_only=(sv_sig and not snp_sig),
                     lead_sv_snp_r2=(round(r2,4) if np.isfinite(r2) else np.nan)))
tab = pd.DataFrame(rows).sort_values(["block","axis"]) if rows else pd.DataFrame()
tab.to_csv(f"{OUTDIR}/nonsnp_peak_dotgrid_lfmm_gifcorr.csv", index=False)
print("wrote", f"{OUTDIR}/nonsnp_peak_dotgrid_lfmm_gifcorr.csv", "|", len(peaks), "peaks x", len(AXES), "axes")''')

nb = new_notebook(); nb["cells"] = C
out = os.path.join(HERE, "nonsnp_peak_dotgrid_lfmm_gifcorr.ipynb")
ep = ExecutePreprocessor(timeout=3600, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": HERE}})
with open(out, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {out}")
