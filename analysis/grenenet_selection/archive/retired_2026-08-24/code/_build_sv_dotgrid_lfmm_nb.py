#!/usr/bin/env python
"""Build raw-LFMM SV-peak x climate-axis dot grid (companion to raw_manhattan_snp_vs_sv_lfmm).

Rows = SV peak blocks (a clq0.9 block whose lead SV record clears per-record Bonferroni
on >=1 axis), gene-labelled. Columns = the 20 climate axes (bio1-19 + pc1). One dot per
(block, axis):
  - dot size  proportional to SV -log10p (lead record in the block at that axis)
  - dot colour = SNP -log10p at the same block (pale = no SNP signal there = genuinely
    SV-only; dark = a SNP tags the same block)
  - BLACK RING = the SV peak is NOT shared by a SNP peak at that axis (SV lead Bonferroni-
    sig, no SNP record in the block Bonferroni-sig) -> the untagged / new SV signal
  - NO ring = a SNP in the same block is also Bonferroni-sig (SNP tags the SV)
Raw (before-WZA) per-record p-values, LFMM K=16 (no GIF), gen9. Single panel (LFMM only).
Self-locating paths (repo moved 2026-07; see memory). No chart title per repo convention.
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor
import os

HERE = os.path.dirname(os.path.abspath(__file__))
C = []
md = lambda s: C.append(new_markdown_cell(s))
co = lambda s: C.append(new_code_cell(s))

md("""# SV peaks x climate axes — raw-LFMM dot grid (SNP-tagging view)

Companion to `raw_manhattan_snp_vs_sv_lfmm.ipynb`, same raw (before-WZA) per-record
LFMM K=16 p-values. **Rows** = SV peak blocks (lead SV record clears per-record
Bonferroni on >=1 axis), gene-labelled. **Columns** = the 20 climate axes.

Per (block, axis): **dot size** = SV -log10p (block lead); **colour = max LD r2 of the
lead SV to any SNP within +/-50 kb** (from the masked `sv_snp_ld_v2` tagging analysis:
dark = the SV is in tight LD with a nearby SNP -> taggable; pale = no nearby SNP tracks
it; **grey = untestable**, anchor too missing). **Black ring** = SV peak here **not
shared by a SNP peak** (SV Bonf-sig, no SNP in the block Bonf-sig).

Two independent "is it a SNP?" axes: the **ring** = does a *climate-significant* SNP sit
in the same block; the **colour** = is the SV in *LD* with any nearby SNP at all
(regardless of that SNP's climate signal). A **pale ringed** dot is the strongest
genuinely-novel-SV candidate: SV-significant, no co-significant SNP, **and** not even
taggable by a nearby SNP. Raw p is uncalibrated (no GIF) - read alongside the per-axis
QQ/lambda in the Manhattan notebook.""")

co('''import os, sys
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), "..", "..")))
import lib
GEA = lib.GEA
WZAIN = f"{GEA}/phase1_replication/results/multiaxis/wza_in"
OUTDIR = f"{GEA}/phase1_replication/results/multiaxis"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
plt.rcParams.update({"figure.dpi": 200})

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

co('''# --- SV -> nearby-SNP LD tagging (best r2 within +/-50 kb), from the tagging analysis ---
# sv_snp_ld_v2/tagging_sv_panel_Chr*.npz : one row per SV anchor, field best_r2 =
# max LD r2 to any SNP within 50 kb (founder-panel genotype space, the same 231-founder
# arch3 panel the GEA decomposes). Also load the indel file so lead records the GEA calls
# "sv" but tagging classifies as indel (size<=50 bp) are still covered. Keyed (chrom,pos,size).
LDDIR = f"{GEA}/sv_snp_ld_v2_maconly"   # MAC floor >=2 kept, missingness gate DROPPED (min-pair-frac 0.0):
                                        # every SV gets an r2 (no untestable/grey), at the cost that r2 from
                                        # few jointly-called founders is noisier. Chosen 2026-07-22 to see all SVs.
R2 = {}          # (chrom, pos, size) -> best_r2 ; POS2 = (chrom,pos) -> max best_r2 fallback
POS2 = {}
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
    return POS2.get((chrom, int(pos)), np.nan)     # fall back to any record at that pos, else untestable
print(f"tagging r2: {len(R2)} anchor records, {len(POS2)} testable positions")

def snp_lead_by_block(axis):
    d = pd.read_csv(f"{WZAIN}/lfmm_snp_gen9_{axis}.csv", usecols=["block", "pval"])
    d = d[d["block"].notna() & (d["block"] != "")]
    bonf = -np.log10(0.05 / len(d))
    return (-np.log10(d["pval"].clip(lower=1e-300))).groupby(d["block"]).max(), bonf

def sv_lead_by_block(axis):
    """lead (max -log10p) SV RECORD per block: keep its chrom/pos/size so we can look up r2."""
    d = pd.read_csv(f"{WZAIN}/lfmm_sv_gen9_{axis}.csv",
                    usecols=["chrom", "pos", "ref_len", "alt_len", "block", "pval"])
    d = d[d["block"].notna() & (d["block"] != "")]
    d["nlp"] = -np.log10(d["pval"].clip(lower=1e-300))
    bonf = -np.log10(0.05 / len(d))
    lead = d.loc[d.groupby("block")["nlp"].idxmax()].copy()
    lead["size"] = (lead["alt_len"] - lead["ref_len"]).abs()
    return lead, bonf

# grid[(block, axis)] = (sv_nlp, sv_sig, snp_sig, r2)  ; r2 = lead SV record's LD tag to a nearby SNP
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
peaks = sorted(recur, key=lambda b: (recur[b], b), reverse=True)   # most-recurrent first (-> bottom row)
print(f"{len(peaks)} SV peak blocks (lead SV Bonferroni-sig on >=1 of 20 axes)")

CAP = None            # show ALL SV peak blocks (no cap), packed dense
if CAP is not None and len(peaks) > CAP:
    print(f"NOTE: {len(peaks)} peaks > {CAP}; showing the {CAP} most recurrent (dropped "
          f"{len(peaks)-CAP} that were SV-sig on the fewest axes).")
    peaks = peaks[:CAP]

def rlabel(b):
    gs = genes_on(b); s = SPANS.loc[b] if b in SPANS.index else None
    pos = f"{s.chrom}:{s.start/1e6:.1f}" if s is not None else b
    return f"{gs[0] if gs else b} {pos}"     # single line, compact

norm = mpl.colors.Normalize(vmin=0, vmax=1.0); cmap = plt.cm.OrRd   # colour = LD r2 (0..1)
NA_COLOR = "0.82"   # grey = SV untestable for tagging (anchor missingness); r2 undefined

def draw(peaks_sub, pitch, dot_base, dot_scale, yfont, ringlw):
    fig, ax = plt.subplots(figsize=(9.0, 0.7 + pitch * len(peaks_sub)))
    n_ring = n_na = 0
    for xi, axis in enumerate(AXES):
        for yi, b in enumerate(peaks_sub):
            g = grid.get((b, axis))
            if g is None: continue
            sv_nlp, sv_sig, snp_sig, r2 = g
            ring = sv_sig and not snp_sig; n_ring += int(ring)
            col = NA_COLOR if not np.isfinite(r2) else cmap(norm(r2)); n_na += int(not np.isfinite(r2))
            ax.scatter(xi, yi, s=dot_base + dot_scale * min(sv_nlp, 10), c=[col],
                       edgecolors=("black" if ring else "none"), linewidths=(ringlw if ring else 0), zorder=3)
    ax.set_xticks(range(len(AXES))); ax.set_xticklabels(AXES, rotation=90, fontsize=5)
    ax.set_yticks(range(len(peaks_sub))); ax.set_yticklabels([rlabel(b) for b in peaks_sub], fontsize=yfont)
    ax.tick_params(length=1.5, pad=1.0)
    ax.set_xlabel("climate axis", fontsize=7); ax.set_ylim(-0.6, len(peaks_sub) - 0.4); ax.margins(x=0.02)
    ax.grid(True, alpha=0.13, lw=0.3, zorder=0)
    ax.annotate(f"raw LFMM \\u00b7 SV \\u2014 size \\u221d SV \\u2212log10p    \\u25cb = SV peak not shared by a SNP peak "
                f"({n_ring} ringed / {len(peaks_sub)} peaks)", (0.0, 1.004), xycoords="axes fraction", fontsize=5.5, va="bottom")
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.03, pad=0.01)
    cb.set_label("max LD r\\u00b2 of the SV to a SNP within 50 kb\\n(dark = tagged by a nearby SNP)", fontsize=6)
    cb.ax.tick_params(labelsize=5, length=2)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    fig.tight_layout(); plt.show()
    print(f"{len(peaks_sub)} peaks: {n_ring} ringed cells; {n_na} grey (untestable)")''')

# ---- Section: all SV peaks (dense) ----
md("""## All SV peaks (dense)

Every SV peak block (lead SV Bonferroni-sig on >=1 axis), most-recurrent at the bottom.""")
co("draw(peaks, pitch=0.125, dot_base=2, dot_scale=7, yfont=3.0, ringlw=0.55)")

# ---- Section: top-45 most recurrent (compact) ----
md("""## Top 45 most-recurrent SV peaks (compact)

Same grid restricted to the 45 SV peaks significant on the most axes — the readable
summary view of the busiest rows above.""")
co("draw(peaks[:45], pitch=0.20, dot_base=3, dot_scale=11, yfont=4.5, ringlw=0.8)")

# ---- export the full underlying table ----
co('''rows = []
for (b, axis), (sv_nlp, sv_sig, snp_sig, r2) in grid.items():
    if b not in peaks: continue
    gs = genes_on(b)
    rows.append(dict(block=b, axis=axis, gene=(gs[0] if gs else ""), sv_nlp=round(sv_nlp,3),
                     sv_bonf_sig=sv_sig, snp_bonf_sig=snp_sig, sv_only=(sv_sig and not snp_sig),
                     lead_sv_snp_r2=(round(r2,4) if np.isfinite(r2) else np.nan)))
tab = pd.DataFrame(rows).sort_values(["block","axis"])
tab.to_csv(f"{OUTDIR}/sv_peak_dotgrid_lfmm.csv", index=False)
print("wrote", f"{OUTDIR}/sv_peak_dotgrid_lfmm.csv", "|", len(peaks), "peaks x", len(AXES), "axes")''')

nb = new_notebook(); nb["cells"] = C
out = os.path.join(HERE, "sv_peak_dotgrid_lfmm.ipynb")
ep = ExecutePreprocessor(timeout=3600, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": HERE}})
with open(out, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {out}")
