#!/usr/bin/env python
"""Build + execute `kmate_benchmarks.ipynb`: the kMate main figure (paper_figures/panel_figure_kmate/).

Inputs are precomputed by `benchmarks/simulations/` (runs.py lists every pool and fit; score.py writes the tables to
`benchmarks/simulations/results/`), so the notebook only reads tables and draws: a rebuild takes seconds.
Panel A: the data sources are drawn here; the method schematic is an Illustrator drawing (`panel_figure_kmate/kMate_Figure1.ai`):
the figure leaves its slot empty at the exact size, to paste it in; `figure_kmate_preview` shows the PNG in the slot,
for viewing only.
Env: basic.  Run: python _build_kmate_benchmarks_nb.py
"""
import os, sys, subprocess
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
NB = f"{HERE}/kmate_benchmarks.ipynb"

md_title = """# kMate benchmarks: the main-figure panels

Accuracy and speed of kMate's production estimator (`--unit chrom`, per-founder normalisation) on simulated pools
built from the 231-founder panel (chromosome 1, 10× unless stated). Data: `benchmarks/simulations/` (`runs.py`: every
pool and fit by figure; `simulate.py`, `fit.py`, `submit.py`; `score.py` → `benchmarks/simulations/results/`). This notebook:
`notebooks/_build_kmate_benchmarks_nb.py` (edit the builder, then rebuild; seconds).

- **Allele frequencies**: every Chr1 panel record passing the project's variant filter (called in at least 50% of
  the 231 founders, minor-allele count at least 2), truth and estimate paired row by row. SNPs, short variants (under 50 bp) and SV alleles (50 bp or more).
- **Founder frequencies and speed**: the pool-size × depth runs (5 seeds each), kMate only. The comparison with
  hapFIRE is supplementary.
- **Simulated genomes**: each founder's Chr1 as the panel VCF records it (consensus on TAIR10), with its missing
  genotypes written as N (`bcftools consensus -M N`, `benchmarks/p231/scripts/05n_build_fastas_p231n.sh`): the panels
  give a founder no k-mer at its missing records, so the reads carry no sequence the panel cannot place. Recombinant
  segments are cut from each founder at its own coordinates (the consensus chain), so every individual keeps its whole
  sequence; reads (wgsim, seeded) are proportional to each individual's length without N."""

code_load = r"""
import os, sys, numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, os.getcwd()); import paper_style as PS; TH = PS.TH
D = "/global/scratch/projects/fc_moilab/tbellg/kmate/benchmarks/simulations/results"
SCHEM = f"{os.getcwd()}/paper_figures/panel_figure_kmate/Copy of kMate_Figure1.png"
SCEN = [("n231_g0", "231 founders"), ("n50_g0", "50 founders"),
        ("n231_g1_self97", "231 founders, 1 generation"), ("n50_g1_self97", "50 founders, 1 generation"),
        ("n50_g3_self97", "50 founders, 3 generations"), ("n50_g3_dom50_self97", "50 founders, 3 generations,\none founder sweeping")]
CLS = [("SNP", "SNPs", TH.GREENS[0]), ("short", "Short variants", TH.GREENS[1]), ("SV", "SV alleles", TH.GREENS[2])]
PTS = {sc: pd.read_csv(f"{D}/af_points_{sc}.tsv.gz", sep="\t") for sc in ["n50_g3_self97"]}
ACC = pd.read_csv(f"{D}/af_accuracy.tsv", sep="\t")
HS = pd.read_csv(f"{D}/pool_size_depth.tsv", sep="\t")       # both tools; the main figure draws kMate's columns
TM = pd.read_csv(f"{D}/timing.tsv", sep="\t")                # speed: 50-founder pools, each run alone on a whole node
print("allele-frequency accuracy (call rate >= 50%, MAC >= 2):")
print(ACC.assign(MAE=ACC.MAE.round(4), RMSE=ACC.RMSE.round(4), R2=ACC.R2.round(4), slope=ACC.slope.round(3)).to_string(index=False))
"""

md_fig = """## Paper figure

`panel_figure_kmate/figure_kmate` (174 × 190 mm): **A** where the data come from (drawn here: long-read assemblies and
short reads → the 231-founder panel; pool-seq → the reads) above a slot for the method schematic (pasted in Illustrator); **B** founder-frequency error (RMSE) against
the number of founders in the pool, at 1×, 5× and 10× (points: seeds; line: median; dotted: 1/N, the average
founder's frequency); **C** true vs estimated allele frequency of SNPs, short variants and SV alleles in a 50-founder
pool after three generations of 97% selfing, 10× (each point a variant, coloured by the local density of variants,
log, relative to the panel's densest bin; key in the middle panel). `figure_kmate_preview` shows the schematic in its slot (raster; not
for the paper). Caption numbers are printed below."""

code_fig = r"""
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path
MM = PS.MM
W_MM, H_MM = 174, 190
W, H = W_MM * MM, H_MM * MM
# A, top: where the data come from, drawn here (sources -> founder panel; pool-seq -> reads). Counts (sample lists of the
# VCFs): pangenome graph = 135 long-read assemblies (82 of GrENE-Net accessions, of which 78 are panel founders: Set-1 and
# Ped-0 were re-genotyped from short reads, the Ct-1 assembly is ambiguous with Stw-0, No-0 is not a founder; data/flag_list.tsv)
# + 153 short-read founders genotyped on that graph = 231 founders.
# A, below: the method schematic, drawn in Illustrator and pasted into an empty slot (its PNG is 1.27:1; ~8.3 pt text at
# 152 mm: the widest that leaves room for the pool-seq link at the page edge).
X0, SW = 10, 152                                       # left edge of A's drawings (= B's y label), schematic width (mm)
K = SW / 142                                           # the source chart was laid out for a 142 mm schematic
SCHEM_PX = (5332, 4212)
SCHEM_SLOT = (X0, 27, SW, SW * SCHEM_PX[1] / SCHEM_PX[0])
def schem_mm(px, py): return SCHEM_SLOT[0] + px * SCHEM_SLOT[2] / SCHEM_PX[0], SCHEM_SLOT[1] + py * SCHEM_SLOT[2] / SCHEM_PX[0]
LONG_C, SHORT_C, POOL_C, PANEL_C = "#7D63AB", "#50926B", "#CB932B", "#8D949D"   # the drawings' purple / green / orange / grey
C_SCEN = "n50_g3_self97"                               # 50 founders, 3 generations of 97% selfing, 10x (pool p231_n50_g3_self97_s42_10x)
def ax_at(fig, x, y, w, h): return fig.add_axes([x * MM / W, 1 - (y + h) * MM / H, w * MM / W, h * MM / H])   # mm from top-left
def letter(x, y, L): PS.letter(fig_, x * MM, y * MM, L)                                                      # mm from top-left

def sources(fig):
    # one overlay axes in mm (x right, y down) for the source chart and the two links into the schematic
    o = fig.add_axes([0, 0, 1, 1], zorder=5); o.set_xlim(0, W_MM); o.set_ylim(H_MM, 0); o.axis("off"); o.patch.set_alpha(0)
    def box(x, y, w, h, txt, col, fill="white"):
        o.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.6", fc=fill, ec=col, lw=PS.LW))
        o.text(x + w / 2, y + h / 2, txt, ha="center", va="center", linespacing=1.15)
    def arrow(p, q, col, lab=None, lab_xy=None, **kw):
        o.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=6, lw=PS.LW, color=col, shrinkA=0, shrinkB=0, **kw))
        if lab: o.text(*lab_xy, lab, ha="center", va="bottom", color=col)
    def X(v): return X0 + K * v                         # x positions on the 142 mm grid, scaled to the schematic
    y0, h = 2, 9
    box(X(1.5), y0, 31 * K, h, "Long-read assemblies\n135 accessions", LONG_C)
    box(X(41), y0, 23 * K, h, "Pangenome\ngraph", LONG_C)
    box(X(84), y0, 23 * K, h, "Short reads\n153 founders", SHORT_C)
    box(X(111), y0, 28 * K, h, "Pool-seq\nevolved populations", POOL_C)
    arrow((X(32.5), y0 + h / 2), (X(41), y0 + h / 2), LONG_C)
    arrow((X(84), y0 + h / 2), (X(64), y0 + h / 2), SHORT_C, "genotyped on\nthe graph", (X(74), y0 + h / 2 - 0.7))
    # founder panel: 78 founders from the assemblies + 153 genotyped from short reads
    py, ph = y0 + h + 4, 9
    box(X(26), py, 53 * K, ph, "Founder panel: 231 founders\n78 long-read + 153 short-read", PANEL_C, fill="#F0F0F0")
    arrow((X(52.5), y0 + h), (X(52.5), py), LONG_C)
    # panel -> founders-pangenome box (its top edge: schematic px row 64)
    arrow((X(52.5), py + ph), (X(52.5), schem_mm(0, 64)[1] + 0.2), PANEL_C)
    # pool-seq -> round the right of the evolved-population box -> the reads
    xr = schem_mm(5262, 0)[0] + 3; xt, yt = schem_mm(5180, 2000); ys = y0 + h / 2; xs = X(139); r = 1.6
    verts = [(xs, ys), (xr - r, ys), (xr, ys), (xr, ys + r), (xr, yt - r), (xr, yt), (xr - r, yt), (xt, yt)]
    codes = [Path.MOVETO, Path.LINETO, Path.CURVE3, Path.CURVE3, Path.LINETO, Path.CURVE3, Path.CURVE3, Path.LINETO]
    o.add_patch(FancyArrowPatch(path=Path(verts, codes), arrowstyle="-|>", mutation_scale=6, lw=PS.LW, color=POOL_C,
                                shrinkA=0, shrinkB=0, joinstyle="round"))

from matplotlib.colors import LogNorm
def density_scatter(a, x, y, nb=100, n_max=150_000, seed=42):
    # the benchmarks' accuracy-panel look: each point coloured by
    # the 2-D histogram count of its bin (viridis, log), densest drawn last, 150k random points, dashed grey 1:1 line
    H, xe, ye = np.histogram2d(x, y, bins=nb, range=[[0, 1], [0, 1]])
    dens = H[np.clip(np.digitize(x, xe) - 1, 0, nb - 1), np.clip(np.digitize(y, ye) - 1, 0, nb - 1)]
    if len(x) > n_max:
        sel = np.random.default_rng(seed).choice(len(x), n_max, replace=False); x, y, dens = x[sel], y[sel], dens[sel]
    o = np.argsort(dens)
    sm = a.scatter(x[o], y[o], c=dens[o] / dens.max(), s=0.8, cmap="viridis", norm=LogNorm(vmin=1 / H.max(), vmax=1),
                   alpha=0.5, edgecolors="none", rasterized=True)
    a.plot([0, 1], [0, 1], ls=(0, (4, 3)), lw=PS.THIN, color="#888888", alpha=0.8, zorder=3)
    return sm

def draw(preview=False):
    global fig_
    fig = fig_ = plt.figure(figsize=(W, H))
    letter(2, 0, "A")
    if preview:
        a = ax_at(fig, *SCHEM_SLOT); a.imshow(plt.imread(SCHEM)); a.axis("off")
    sources(fig)
    right = X0 + SW + 1                                 # B + C span the drawings' width
    s, gap = 26, 3
    y1 = SCHEM_SLOT[1] + SCHEM_SLOT[3] + 7              # top of B and C's axes (same height, x labels on one line)
    # ---- B: founder-frequency RMSE vs pool size --------------------------------------------------------------
    a = ax_at(fig, X0 + 13, y1, 42, s); letter(2, y1 - 5, "B")
    Ns = sorted(HS.N.unique()); xpos = {n: i for i, n in enumerate(Ns)}
    DEP = [(1, "1×", TH.GREENS[0]), (5, "5×", TH.GREENS[1]), (10, "10×", TH.GREENS[2])]
    for k, (dp, dl, col) in enumerate(DEP):
        d = HS[(HS.depth == dp)]
        x = d.N.map(xpos) + (k - 1) * 0.18
        a.scatter(x, d.kmate_h_RMSE, s=5, color=col, alpha=0.7, linewidths=0, zorder=3)
        m = d.groupby("N").kmate_h_RMSE.median()
        a.plot([xpos[n] + (k - 1) * 0.18 for n in m.index], m.values, color=col, lw=PS.LW, label=dl, zorder=4)
    a.plot(range(len(Ns)), [1 / n for n in Ns], color=TH.MUTED, lw=PS.REF_LW, ls=":", zorder=2)
    a.annotate("1/N", (1, 1 / Ns[1]), xytext=(1, -7), textcoords="offset points", ha="left", va="top", color=TH.MUTED)
    a.set_yscale("log"); a.set_xticks(range(len(Ns))); a.set_xticklabels([str(n) for n in Ns])
    a.set_xlabel("Founders in the pool"); a.set_ylabel("Founder frequency\nerror (RMSE)")
    a.set_xlim(-0.5, len(Ns) - 0.5); a.set_ylim(2e-4, 1.6)
    a.legend(title="Depth", loc="upper right", frameon=False, handlelength=1.2)
    # ---- C: truth vs estimate, one scenario, columns = SNPs / short variants / SV alleles ---------------------
    x0 = right - 3 * s - 2 * gap
    letter(x0 - 11, y1 - 5, "C")
    for c, (cl, clab, _) in enumerate(CLS):
        a = ax_at(fig, x0 + c * (s + gap), y1, s, s)
        d = PTS[C_SCEN][PTS[C_SCEN].cls == cl]
        sm = density_scatter(a, d.truth_af.to_numpy(), d.est_af.to_numpy())
        a.set_xlim(-0.03, 1.03); a.set_ylim(-0.03, 1.03); a.set_aspect("equal")
        tk = [0, 0.5, 1]; a.set_xticks(tk); a.set_yticks(tk)
        a.set_xticklabels(["0", "0.5", "1"]); a.set_yticklabels(["0", "0.5", "1"]); a.tick_params(length=0)
        a.annotate(clab, (0.5, 1.03), xycoords="axes fraction", ha="center", va="bottom")
        if c == 0: a.annotate("50 founders,\n3 generations", (0.04, 0.96), xycoords="axes fraction", ha="left", va="top")
        a.annotate(f"RMSE\n{ACC[(ACC.scenario == C_SCEN) & (ACC.cls == cl)].RMSE.iloc[0]:.4f}", (0.96, 0.04),
                   xycoords="axes fraction", ha="right", va="bottom")
        if c == 1: a.set_xlabel("True allele frequency")
        if c == 0: a.set_ylabel("Estimated allele\nfrequency")
        else: a.set_yticklabels([])
        if c == 1:                                      # colour key in the empty top-left corner of the middle panel
            cax = a.inset_axes([0.05, 0.76, 0.56, 0.045])
            cb = fig.colorbar(sm, cax=cax, orientation="horizontal"); cb.outline.set_visible(False)
            cb.set_ticks([sm.norm.vmin, sm.norm.vmax]); cb.set_ticklabels(["Low", "High"])
            cax.tick_params(which="both", length=0, pad=1); cax.xaxis.set_minor_locator(plt.NullLocator())
            lo, hi = cax.get_xticklabels(); lo.set_ha("left"); hi.set_ha("right")
            cax.annotate("Variant density", (0, 1.8), xycoords="axes fraction", ha="left", va="bottom")
    return fig

with PS.style():
    fig = draw(); PS.save(fig, "panel_figure_kmate", "figure_kmate", tight=False); plt.show()
    fig = draw(preview=True); PS.save(fig, "panel_figure_kmate", "figure_kmate_preview", tight=False); plt.close(fig)

# caption numbers
print("C: MAE, RMSE, R2 per scenario x class: see the table above (scenario n50_g3_self97)")
h = HS.groupby(["N", "depth"]).agg(rmse=("kmate_h_RMSE", "median"), found=("kmate_n_found", "median"), R2=("kmate_h_R2", "median"))
print("\nB: founder-frequency RMSE (median of 5 seeds), founders recovered, R2(h):"); print(h.round(4).to_string())
t = TM.groupby("depth")[["kmate_elapsed_s", "kmate_cpu_s", "kmate_max_rss_mb"]].median()
print("\nSpeed (caption / supplementary): N = 50, alone on a whole node, median wall-clock s, CPU s, peak memory MB:"); print(t.round(0).to_string())
"""


md_supp = """# Supplementary figures (exported)

Each to its own folder under `paper_figures/`, full page width, drawn at print size. Inputs:
`benchmarks/simulations/score.py` → `benchmarks/simulations/results/`. Every
allele-frequency score uses the project variant filter (call rate at least 50%, minor-allele count at least 2)."""

md_s_af = """## Supplementary: allele-frequency accuracy in every scenario

`kmate_af_accuracy/af_accuracy_scenarios`: true vs estimated allele frequency, rows = SNPs / short variants / SV
alleles, columns = the eight simulated pools (1, 2, 231 or 50 founders; 0, 1 or 3 generations of 97% selfing; one
founder sweeping to 50%), 10×, Chr1. In the 1- and 2-founder pools every true frequency is 0, 0.5 or 1, so the panels
show the estimator's spread around those values directly. Points coloured by the local density of variants (log, relative to each panel's densest
bin), dashed 1:1 line; RMSE in each panel (MAE / R²: the table printed under the main figure's data cell)."""

code_s_af = r"""
SCEN8 = [("n1_g0", "1 founder"), ("n2_g0", "2 founders")] + SCEN
PTS6 = {sc: pd.read_csv(f"{D}/af_points_{sc}.tsv.gz", sep="\t") for sc, _ in SCEN8}
def axes_mm(fig, x, y, w, h):
    Wf, Hf = fig.get_size_inches(); return fig.add_axes([x * MM / Wf, 1 - (y + h) * MM / Hf, w * MM / Wf, h * MM / Hf])
HEAD = {"n1_g0": "1 founder", "n2_g0": "2 founders", "n231_g0": "231\nfounders", "n50_g0": "50\nfounders",
        "n231_g1_self97": "231 founders,\n1 gen.", "n50_g1_self97": "50 founders,\n1 gen.", "n50_g3_self97": "50 founders,\n3 gen.",
        "n50_g3_dom50_self97": "50 founders,\n3 gen., one\nsweeping"}
with PS.style():
    s_, g_, x0_, y0_ = 18.2, 1.7, 13, 12                  # 13 + 8 x 18.2 + 7 x 1.7 = 170.5 mm of 174
    fig = plt.figure(figsize=(PS.FULL_W, (y0_ + 3 * s_ + 2 * g_ + 26) * MM))
    for r, (cl, clab, _) in enumerate(CLS):
        for c, (sc, _) in enumerate(SCEN8):
            a = axes_mm(fig, x0_ + c * (s_ + g_), y0_ + r * (s_ + g_), s_, s_)
            d = PTS6[sc][PTS6[sc].cls == cl]
            sm = density_scatter(a, d.truth_af.to_numpy(), d.est_af.to_numpy())
            a.set_xlim(-0.03, 1.03); a.set_ylim(-0.03, 1.03); a.set_aspect("equal"); a.tick_params(length=0)
            tk = [0, 0.2, 0.4, 0.6, 0.8, 1]; a.set_xticks(tk); a.set_yticks(tk)
            a.set_xticklabels(["0", "", "", "", "", "1"] if r == 2 else []); a.set_yticklabels(["0", "", "", "", "", "1"] if c == 0 else [])
            if r == 0: a.annotate(HEAD[sc], (0.5, 1.04), xycoords="axes fraction", ha="center", va="bottom")
            if c == 0: a.annotate(clab, (0.05, 0.95), xycoords="axes fraction", ha="left", va="top")
            rm = ACC[(ACC.scenario == sc) & (ACC.cls == cl)].RMSE.iloc[0]
            a.annotate(f"RMSE\n{rm:.4f}", (0.96, 0.04), xycoords="axes fraction", ha="right", va="bottom")   # two lines: 8 narrow panels
    fig.text((x0_ - 9) * MM / PS.FULL_W, 1 - (y0_ + 1.5 * s_ + g_) * MM / fig.get_size_inches()[1], "Estimated allele frequency",
             rotation=90, ha="center", va="center")
    fig.text((x0_ + 4 * (s_ + g_) - g_ / 2) * MM / PS.FULL_W, 1 - (y0_ + 3 * s_ + 2 * g_ + 6) * MM / fig.get_size_inches()[1],
             "True allele frequency", ha="center", va="top")
    cax = axes_mm(fig, x0_ + 3 * (s_ + g_), y0_ + 3 * s_ + 2 * g_ + 15, 50, 1.6)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal"); cb.outline.set_visible(False)
    cb.set_ticks([sm.norm.vmin, sm.norm.vmax]); cb.set_ticklabels(["Low", "High"])
    cax.tick_params(which="both", length=0, pad=1.5); cax.xaxis.set_minor_locator(plt.NullLocator())
    cax.annotate("Local density of variants (log)", (-0.03, 0.5), xycoords="axes fraction", ha="right", va="center")
    PS.save(fig, "kmate_af_accuracy", "af_accuracy_scenarios", tight=False); plt.show()
"""

md_s_hf = """## Supplementary: kMate and hapFIRE

`kmate_vs_hapfire/kmate_vs_hapfire`: the same simulated pools (1 to 231 founders, 1×, 5×, 10×, 5 seeds, Chr1) run
through kMate (231-founder panel) and hapFIRE (its own GrENE-Net SNP panel); each tool on its own panel, the only fair
basis (a shared panel has to impute missing genotypes as reference for hapFIRE). Columns = depth. **A** founder-frequency
error (RMSE; dotted 1/N). **B** allele-frequency error (RMSE; SNPs only, as for every hapFIRE comparison; the
call-rate-90% score is in the printed table, indistinguishable from this one) on each tool's own SNPs passing the project filter (kMate: call rate at least 50% and minor-allele count at
least 2; hapFIRE's panel is fully imputed: minor-allele count at least 2). **C, D** time per pool, 50 founders, against
depth: wall-clock and total CPU, each run alone on a whole node of one type (savio3, 40 cores; kMate 8 threads; bwa 8
threads; hapFIRE single-threaded in this build, alignment included). Points:
seeds; lines: medians. In A, errors below 1e-5 (hapFIRE in the one-founder pool, about 1e-7) are drawn at 1e-5."""

code_s_hf = r"""
HA = HH = HS                                         # pool_size_depth.tsv: both tools
TOOLS = [("kmate", "kMate", TH.GREENS[2]), ("hapfire", "hapFIRE", "#8A8A8A")]
FLOOR = 1e-5                                         # log panels: values below (hapFIRE's exact 1-founder fits) drawn at the floor
def tool_lines(a, df, col_fmt, Ns, logy=False):
    xp = {n: i for i, n in enumerate(Ns)}
    for k, (t, lab, colr) in enumerate(TOOLS):
        col = col_fmt.format(t); d = df[df.N.isin(Ns) & df[col].notna()].copy()
        if logy: d[col] = d[col].clip(lower=FLOOR)
        off = (k - 0.5) * 0.22
        a.scatter(d.N.map(xp) + off, d[col], s=4, color=colr, alpha=0.75, linewidths=0, zorder=3)
        m = d.groupby("N")[col].median()
        a.plot([xp[n] + off for n in m.index], m.values, color=colr, lw=PS.LW, label=lab, zorder=4)
    a.set_xticks(range(len(Ns))); a.set_xticklabels([str(n) for n in Ns]); a.set_xlim(-0.5, len(Ns) - 0.5)
    if logy: a.set_yscale("log")
with PS.style():
    w_, h_, gx, gy, x0_, y0_ = 44, 30, 9, 13, 19, 8         # 19 + 3 x 44 + 2 x 9 = 169 mm of 174
    Hfig = y0_ + 2 * h_ + gy + 22 + h_ + 12
    fig = plt.figure(figsize=(PS.FULL_W, Hfig * MM))
    Nall = sorted(HH.N.unique()); N150 = [n for n in Nall if n < 231]
    rows_ = [("A", "{}_h_RMSE", Nall, True, "Founder frequency error\n(RMSE)", HH),
             ("B", "{}_af_RMSE", Nall, False, "Allele-frequency error\n(RMSE)", HA)]
    for r, (L, colf, Ns, lg, ylab, df) in enumerate(rows_):
        for c, dp in enumerate([1, 5, 10]):
            a = axes_mm(fig, x0_ + c * (w_ + gx), y0_ + r * (h_ + gy), w_, h_)
            tool_lines(a, df[df.depth == dp], colf, Ns, logy=lg)
            if L == "A":
                a.plot(range(len(Ns)), [1 / n for n in Ns], color=TH.MUTED, lw=PS.REF_LW, ls=":", zorder=2)
                if c == 0: a.annotate("1/N", (1, 1 / Ns[1]), xytext=(4, 2), textcoords="offset points", color=TH.MUTED, va="bottom")
            if r == 0: a.annotate(f"{dp}×", (0.5, 1.06), xycoords="axes fraction", ha="center", va="bottom")
            if c == 0: a.set_ylabel(ylab)
            if r == len(rows_) - 1: a.set_xlabel("Founders in the pool")
            if r == 0 and c == 2: a.legend(loc="upper right", frameon=False, handlelength=1.2)
        PS.letter(fig, 2 * MM, (y0_ + r * (h_ + gy) - 2) * MM, L)
    # share y within each row
    for r in range(len(rows_)):
        ax_r = fig.axes[3 * r: 3 * r + 3]; lo = min(a.get_ylim()[0] for a in ax_r); hi = max(a.get_ylim()[1] for a in ax_r)
        for a in ax_r: a.set_ylim((FLOOR / 1.6, 2) if rows_[r][3] else (lo, hi))
    yS = y0_ + 2 * h_ + gy + 22
    d50 = TM                                              # timing runs: each alone on a whole node
    for k, (L, kcol, hcol, ylab) in enumerate([("C", "kmate_elapsed_s", "hapfire_total_elapsed_s", "Wall-clock time per pool (min)"),
                                               ("D", "kmate_cpu_s", "hapfire_total_cpu_s", "CPU time per pool (min)")]):
        a = axes_mm(fig, x0_ + k * 80, yS, 62, h_); PS.letter(fig, (2 + k * 80) * MM, (yS - 4) * MM, L)
        dd = d50[d50[kcol].notna() & d50[hcol].notna()]; dps = sorted(dd.depth.unique()); xp = {v: i for i, v in enumerate(dps)}
        for j, (col, lab, colr) in enumerate([(kcol, "kMate", TH.GREENS[2]), (hcol, "hapFIRE", "#8A8A8A")]):
            off = (j - 0.5) * 0.2
            a.scatter(dd.depth.map(xp) + off, dd[col] / 60, s=4, color=colr, alpha=0.75, linewidths=0, zorder=3)
            m = dd.groupby("depth")[col].median() / 60
            a.plot([xp[v] + off for v in m.index], m.values, color=colr, lw=PS.LW, label=lab, zorder=4)
        a.set_yscale("log"); a.set_xticks(range(len(dps))); a.set_xticklabels([f"{int(v)}×" for v in dps])
        a.yaxis.set_major_locator(plt.FixedLocator([1, 2, 5, 10, 20, 50])); a.yaxis.set_minor_locator(plt.NullLocator())
        a.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
        a.set_xlabel("Sequencing depth"); a.set_ylabel(ylab)
        if k == 0: a.legend(loc="center left", frameon=False, handlelength=1.2)
    PS.save(fig, "kmate_vs_hapfire", "kmate_vs_hapfire", tight=False); plt.show()
print("allele-frequency accuracy, median over seeds (project filter):")
print(HA.groupby(["N", "depth"])[["kmate_af_RMSE", "kmate_af_RMSE_cr90", "hapfire_af_RMSE", "kmate_n_snps", "kmate_n_snps_cr90", "hapfire_n_snps"]].median().round(4).to_string())
print("\nfounder frequencies, median over seeds:")
print(HH.groupby(["N", "depth"])[["kmate_h_R2", "hapfire_h_R2", "kmate_h_RMSE", "hapfire_h_RMSE", "kmate_n_found", "hapfire_n_found",
                                  ]].median().round(4).to_string())
print("\ntime per pool, 50 founders, alone on a whole node, median over seeds (s):")
print(TM.groupby("depth")[["kmate_elapsed_s", "hapfire_total_elapsed_s", "kmate_cpu_s", "hapfire_total_cpu_s"]].median().round(0).to_string())
"""

md_s_ft = """## Supplementary: the mixed-panel imbalance and how kMate removes it

`kmate_founder_type/founder_type`: estimated frequency of every founder (crosses) against its true frequency (grey),
founders from long-read assemblies (dark, the 78 of the 231-founder panel) and from short-read genotyping (light, 153), in
a pool of all 231 founders at equal frequency and in two 50-founder pools weighted to one founder type (40 long-read + 10
short-read; 45 short-read + 5 long-read), 10×, Chr1, seed 42. Rows add kMate's two corrections in turn. **A** all panel
k-mers, every founder normalized by the same total (the multinomial update that ignores each founder's k-mer content).
**B** private k-mers removed (carried by one founder only; k-mers carried by every founder also dropped), same
normalization. **C** private k-mers removed and each founder normalized by its own k-mer content (production). In the
50-founder pools founders are ordered by true frequency. RMSE per panel, the mean error per founder type and the number
of founders estimated below 0.001: printed below (caption). Data: figure `founder_type` in `benchmarks/simulations/runs.py`
→ `founder_type.tsv`. Founder type follows the panel: long-read = the 78 samples of the cactus VCF."""

code_s_ft = r"""
FT = pd.read_csv(f"{D}/founder_type.tsv", sep="\t")
REG = [("n231_g0", "231 founders, equal frequency", False), ("n50_longheavy", "50 founders: 40 long-read, 10 short-read", True),
       ("n50_shortheavy", "50 founders: 45 short-read, 5 long-read", True)]
ARM = [("unfilt_global", "All k-mers,\nshared normalization"), ("filt_global", "Private k-mers removed,\nshared normalization"),
       ("filt", "Private k-mers removed,\nper-founder normalization")]
C_LONG, C_SHORT = TH.GREENS[2], TH.GREENS[0]
with PS.style():
    w_, h_, gx, gy, x0_, y0_ = 44, 30, 7, 9, 22, 8          # 22 + 3 x 44 + 2 x 7 = 168 mm of 174
    fig = plt.figure(figsize=(PS.FULL_W, (y0_ + len(ARM) * h_ + (len(ARM) - 1) * gy + 22) * MM)); rows_ = []
    for r, (arm, alab) in enumerate(ARM):
        for c, (reg, rlab, by_truth) in enumerate(REG):
            a = axes_mm(fig, x0_ + c * (w_ + gx), y0_ + r * (h_ + gy), w_, h_)
            d = FT[(FT.regime == reg) & (FT.arm == arm)].reset_index(drop=True)
            o = np.argsort(d.truth.to_numpy(), kind="stable") if by_truth else np.arange(len(d)); d = d.iloc[o].reset_index(drop=True)
            x = np.arange(len(d))
            a.scatter(x, d.truth, s=3, color="#8A8A8A", linewidths=0, zorder=1, label="True frequency")
            for lr, colr, lab in ((False, C_SHORT, "Short-read genotyping founders"), (True, C_LONG, "Long-read assembly founders")):
                m = d.long_read == lr
                a.scatter(x[m], d.h[m], marker="x", s=6, linewidths=0.5, color=colr, zorder=2 + lr, label=lab)
            if c == 0: PS.panel(a, "ABC"[r])
            if r == 0: a.annotate(rlab, (0.5, 1.05), xycoords="axes fraction", ha="center", va="bottom")
            if c == 0:
                a.set_ylabel("Founder frequency")
                a.annotate(alab, (-0.42, 0.5), xycoords="axes fraction", ha="center", va="center", rotation=90)
            if r == len(ARM) - 1: a.set_xlabel("Founders (ordered by true frequency)" if by_truth else "Founders")
            else: a.set_xticklabels([])
            a.set_ylim(-0.003, 0.07 if c else 0.05)
    for arm, alab in ARM + [("unfilt", "All k-mers, per-founder normalization (not drawn)")]:   # every arm, for the caption
        for reg, rlab, _ in REG:
            d = FT[(FT.regime == reg) & (FT.arm == arm)]; e = d.h - d.truth
            rows_.append(dict(row=alab.replace("\n", " "), pool=rlab, RMSE=np.sqrt((e ** 2).mean()),
                              mean_err_long=e[d.long_read].mean(), mean_err_short=e[~d.long_read].mean(),
                              n_below_0_001=int((d.h < 1e-3).sum())))
    h1, l1 = fig.axes[0].get_legend_handles_labels()
    fig.legend(h1[::-1], l1[::-1], loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.2)
    PS.save(fig, "kmate_founder_type", "founder_type", tight=False); plt.show()
print(pd.DataFrame(rows_).round(5).to_string(index=False))
"""


md_s_lr = """## Supplementary: a panel of long-read founders only

`kmate_longread_panel/longread_panel`: pools drawn only from the 78 long-read assembly founders that are in every panel
(the 80 of the long-read panel minus 6910 and 7275, absent from the 231 panel and the GrENE-Net VCF), 1 to 78 founders at
equal frequency, generation 0, 1×, 5×, 10×, 5 seeds, Chr1; reads simulated from each founder's genome as the long-read
panel VCF records it (consensus on TAIR10, `bcftools consensus -M N`: a founder's missing genotypes are written as N, as the
production k-mer panel treats them, so the reads carry no sequence the panel cannot place; for hapFIRE, from its fully
imputed GrENE-Net SNP genomes, same founder draw). Arms: kMate with the production 231-founder panel on the same
reads; kMate with a 78-founder panel (`benchmarks/p80/scripts/build_p78_panels.py`: records and k-mers no longer
variable among the 78 dropped) with and without the private k-mer filter (recomputed for 78 founders); hapFIRE with the
GrENE-Net VCF cut to the 78 (sites no longer polymorphic dropped). The figure shows kMate's best setting
(long-read panel, private k-mer filter) against hapFIRE; every arm is in the tables below. Columns = depth. **A**
founder-frequency error (RMSE over the 78 founders; below 1e-5 drawn at 1e-5). **B** allele-frequency error (RMSE) on
SNPs only (as for every hapFIRE comparison) with minor-allele count at least 2 among the 78 founders, kMate: called in at
least 50% of the 78 (at least 90%: in the tables below, indistinguishable; hapFIRE's panel is fully imputed). Data: figure
`longread_grid` in `benchmarks/simulations/runs.py` (kMate and hapFIRE) → `longread_grid.tsv`."""

code_s_lr = r"""
LR = pd.read_csv(f"{D}/longread_grid.tsv", sep="\t")
ARMS_ALL = [("p231", "kMate, 231-founder panel", TH.GREENS[0], "-"), ("p78_unfilt", "kMate, long-read panel, no filter", TH.GREENS[1], "-"),
            ("p78_filt", "kMate, long-read panel, private k-mer filter", TH.GREENS[2], "-"), ("hapfire", "hapFIRE, long-read founders", "#8A8A8A", "-")]
ARMS = [("p78_filt", "kMate", TH.GREENS[2], "-"), ("hapfire", "hapFIRE", "#8A8A8A", "-")]   # plotted: kMate's best setting vs hapFIRE
NL = sorted(LR.N.unique()); xpl = {n: i for i, n in enumerate(NL)}
def arm_lines(a, df, col, logy=False):
    for k, (arm, lab, colr, ls) in enumerate(ARMS):
        d = df[(df.arm == arm) & df[col].notna()].copy()
        if logy: d[col] = d[col].clip(lower=FLOOR)
        off = (k - 0.5) * 0.22
        a.scatter(d.N.map(xpl) + off, d[col], s=4, color=colr, alpha=0.75, linewidths=0, zorder=3)
        m = d.groupby("N")[col].median()
        a.plot([xpl[n] + off for n in m.index], m.values, color=colr, lw=PS.LW, ls=ls, label=lab, zorder=4)
    a.set_xticks(range(len(NL))); a.set_xticklabels([str(n) for n in NL]); a.set_xlim(-0.5, len(NL) - 0.5)
    if logy: a.set_yscale("log")
with PS.style():
    w_, h_, gx, gy, x0_, y0_ = 44, 30, 9, 13, 19, 8
    fig = plt.figure(figsize=(PS.FULL_W, (y0_ + 2 * h_ + gy + 21) * MM))
    for r, (L, col, lg, ylab) in enumerate([("A", "h_rmse78", True, "Founder frequency error\n(RMSE)"),
                                            ("B", "af_rmse", False, "Allele-frequency error\n(RMSE)")]):
        axr = []
        for c, dp in enumerate([1, 5, 10]):
            a = axes_mm(fig, x0_ + c * (w_ + gx), y0_ + r * (h_ + gy), w_, h_); axr.append(a)
            arm_lines(a, LR[LR.depth == dp], col, logy=lg)
            if r == 0: a.annotate(f"{dp}×", (0.5, 1.06), xycoords="axes fraction", ha="center", va="bottom")
            if c == 0: a.set_ylabel(ylab)
            if r == 1: a.set_xlabel("Founders in the pool")
        lo = min(a.get_ylim()[0] for a in axr); hi = max(a.get_ylim()[1] for a in axr)
        for a in axr: a.set_ylim((FLOOR / 1.6, 0.05) if lg else (0, hi))
        PS.letter(fig, 2 * MM, (y0_ + r * (h_ + gy) - 2) * MM, L)
    hl, ll = fig.axes[5].get_legend_handles_labels()                   # below the panels: the 10x panel has no free corner
    fig.legend(hl, ll, loc="lower center", ncol=3, frameon=False, handlelength=2.0, bbox_to_anchor=(0.5, 0.0))
    PS.save(fig, "kmate_longread_panel", "longread_panel", tight=False); plt.show()
lab = {a: l for a, l, _, _ in ARMS_ALL}
for col, nm in (("h_rmse78", "founder-frequency RMSE (78 founders)"), ("h_misassigned", "frequency put on the wrong founders (0.5 x sum |h - truth|)"),
                ("af_rmse", "allele-frequency RMSE, SNPs, call rate >= 50%"), ("af_rmse_cr90", "allele-frequency RMSE, SNPs, call rate >= 90%"),
                ("af_rmse_all", "allele-frequency RMSE, SNPs + short variants + SVs, call rate >= 50% (kMate only)")):
    print(f"\n{nm}, 10x, median of 5 seeds:")
    print(LR[LR.depth == 10].pivot_table(index="N", columns="arm", values=col, aggfunc="median").reindex(columns=[a for a, *_ in ARMS_ALL]).dropna(axis=1, how="all").rename(columns=lab).round(4).to_string())
print("\nfrequency on short-read founders, kMate 231 panel (median):"); print(LR[LR.arm == "p231"].pivot_table(index="N", columns="depth", values="h_on_short", aggfunc="median").round(3).to_string())
print("\none-founder pools, frequency on the true founder (median over depths and seeds):"); print(LR[LR.N == 1].groupby("arm").h_true_founder.median().rename(index=lab).round(3).to_string())
print("\ntime per pool, 50 founders (min, median):"); print((LR[LR.N == 50].groupby(["arm", "depth"])[["elapsed_s", "cpu_s"]].median() / 60).round(1).to_string())
"""

md_ub = """# Estimation units: whole chromosome, windows and LD blocks (exploratory, not exported)

kMate fitted once per chromosome (`--unit chrom`), in fixed windows (`--unit bp --window-bp N`) or in r² ≥ 0.1 LD
blocks of the panel (`--unit ld`), on pools of 50 founders (Chr1, 10×) after 1, 3 or 10 generations of outcrossing or of
97% selfing, 3 seeds each. **Long-read panel**: pools drawn from the 78 long-read assembly founders, fitted with the
78-founder panel without the private k-mer filter (every k-mer kept). **Mixed panel**: pools drawn from the 231 founders,
fitted with the production 231-founder panel. Variant filter as in the main figure (called in at least 50% of the
panel's founders, minor-allele count at least 2). A window holding fewer than 200 observed k-mers gives no estimate; errors
are over the variants estimated, and the share estimated is drawn separately. Scored whole-chromosome fits on exactly
the variants each window estimated give the same error as on all variants, so skipped windows do not flatter the window
units. The LD partition is coarse (5 blocks on Chr1 for the long-read panel, median 5.7 Mb; 18 for the mixed panel).

Runs: figures `units_p78` and `units_p231` in `benchmarks/simulations/runs.py` → `units.tsv`. Figures are shown here
and written to `notebooks/plots/kmate_unitbench_*.png`."""

code_ub_load = r"""
U = pd.read_csv(f"{D}/units.tsv", sep="\t")
UB, UBM = U[U.founders == "p78"], U[U.founders == "p231"]           # long-read panel, mixed panel
def wkb(u): return int(u[2:-1]) if u.startswith("bp") else None
WINS = sorted({u for u in UB.unit if u.startswith("bp")}, key=wkb)
UNITS_UB = ["chrom"] + WINS + ["ld"]
WCOL = dict(zip(WINS, plt.cm.viridis(np.linspace(0.05, 0.85, len(WINS)))))
ULAB = {"chrom": "Whole chromosome", "ld": "LD partition (r2 0.1)", **{u: f"{wkb(u)} kb windows" for u in WINS}}
UCOL = {"chrom": "black", "ld": "#8A8A8A", **WCOL}
ULS = {"ld": (0, (3, 2))}
MAT = [("outcrossing", "Outcrossing"), ("selfing 97%", "97% selfing")]
GENS = [1, 3, 10]; xg = {g: i for i, g in enumerate(GENS)}
PLOTS = f"{os.getcwd()}/plots"
def mean_tab(df, val, cls="SNP"):
    return (df[df.cls == cls].groupby(["mating", "generations", "unit"])[val].mean().unstack("unit")
            .reindex(columns=[u for u in UNITS_UB if u in set(df.unit)]))
print("\nSNP allele-frequency RMSE, mean of 3 seeds:"); print(mean_tab(UB, "RMSE").round(4).to_string())
print("\nshare of SNPs estimated:"); print(mean_tab(UB, "frac_estimated").round(3).to_string())
"""

md_ub_gen = """## Error against generations, per unit and variant class (long-read panel)

Rows: SNPs, short variants, SV alleles. Columns: outcrossing, 97% selfing. Allele-frequency RMSE (log) after 1, 3 and
10 generations; points are seeds, lines the mean of 3 seeds."""

code_ub_gen = r"""
with PS.style():
    fig, axs = plt.subplots(3, 2, figsize=(PS.FULL_W, 150 * MM), sharex=True, sharey="row")
    for r, (cl, clab, _) in enumerate(CLS):
        for c, (mt, mlab) in enumerate(MAT):
            a = axs[r, c]
            for k, u in enumerate(UNITS_UB):
                d = UB[(UB.cls == cl) & (UB.mating == mt) & (UB.unit == u)]
                if d.empty: continue
                a.scatter(d.generations.map(xg), d.RMSE, s=4, color=UCOL[u], alpha=0.6, linewidths=0, zorder=3)
                m = d.groupby("generations").RMSE.mean()
                a.plot([xg[g] for g in m.index], m.values, color=UCOL[u], ls=ULS.get(u, "-"), lw=PS.LW, label=ULAB[u], zorder=4)
            a.set_yscale("log"); a.set_xticks(range(len(GENS))); a.set_xticklabels([str(g) for g in GENS]); a.set_xlim(-0.15, len(GENS) - 0.85)
            a.annotate(f"{clab}, {mlab.lower()}", (0.03, 0.95), xycoords="axes fraction", ha="left", va="top")
            if c == 0: a.set_ylabel("Allele-frequency error\n(RMSE)")
            if r == 2: a.set_xlabel("Generations")
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, frameon=False, handlelength=2.0, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.17, hspace=0.12, wspace=0.08)
    fig.savefig(f"{PLOTS}/kmate_unitbench_generations.png", dpi=200); plt.show()
"""

md_ub_sweep = """## Window size (long-read panel, SNPs)

**Top**: allele-frequency RMSE against window size (log), one line per number of generations; the dashed line of the same
colour is the whole-chromosome fit of those pools, the open circle each line's best window. **Bottom**: the share of SNPs
estimated (windows with fewer than 200 observed k-mers give none). Points are seeds, lines the mean of 3 seeds."""

code_ub_sweep = r"""
GCOL = dict(zip(GENS, TH.GREENS[:3]))
best = []
with PS.style():
    fig, axs = plt.subplots(2, 2, figsize=(PS.FULL_W, 110 * MM), sharex=True, sharey="row", gridspec_kw=dict(height_ratios=[2, 1]))
    xs = [wkb(u) for u in WINS]
    for c, (mt, mlab) in enumerate(MAT):
        S = UB[(UB.cls == "SNP") & (UB.mating == mt)]
        for g in GENS:
            d = S[(S.generations == g) & S.unit.isin(WINS)].assign(kb=lambda x: x.unit.map(wkb))
            m = d.groupby("kb")[["RMSE", "frac_estimated"]].mean().reindex(xs)
            axs[0, c].scatter(d.kb, d.RMSE, s=4, color=GCOL[g], alpha=0.6, linewidths=0, zorder=3)
            axs[0, c].plot(m.index, m.RMSE, color=GCOL[g], lw=PS.LW, marker="o", ms=2, label=f"{g} generation{'s' * (g > 1)}", zorder=4)
            ch = S[(S.generations == g) & (S.unit == "chrom")].RMSE.mean()
            axs[0, c].axhline(ch, color=GCOL[g], ls=(0, (4, 3)), lw=PS.REF_LW, zorder=2)
            kb_best = m.RMSE.idxmin()
            axs[0, c].scatter([kb_best], [m.RMSE.min()], s=40, facecolors="none", edgecolors=GCOL[g], linewidths=PS.THIN, zorder=5)
            axs[1, c].scatter(d.kb, d.frac_estimated, s=4, color=GCOL[g], alpha=0.6, linewidths=0, zorder=3)
            axs[1, c].plot(m.index, m.frac_estimated, color=GCOL[g], lw=PS.LW, zorder=4)
            best.append(dict(mating=mt, generations=g, best_window_kb=kb_best, best_window_RMSE=m.RMSE.min(), chrom_RMSE=ch,
                             frac_estimated_at_best=m.frac_estimated.loc[kb_best]))
        axs[0, c].annotate(mlab, (0.03, 0.95), xycoords="axes fraction", ha="left", va="top")
        axs[1, c].set_xlabel("Window size (kb)")
    for a in axs.ravel():
        a.set_xscale("log"); a.set_xticks(xs); a.set_xticklabels([str(x) for x in xs]); a.xaxis.set_minor_locator(plt.NullLocator())
    axs[0, 0].set_yscale("log"); axs[0, 0].set_ylabel("Allele-frequency error\n(RMSE)"); axs[1, 0].set_ylabel("Share of SNPs\nestimated")
    axs[1, 0].set_ylim(0, 1.02)
    axs[0, 1].legend(loc="upper right", frameon=False, handlelength=1.5)
    fig.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.1, hspace=0.08, wspace=0.08)
    fig.savefig(f"{PLOTS}/kmate_unitbench_window_sweep.png", dpi=200); plt.show()
print("best window per condition (SNP RMSE, mean of 3 seeds) against the whole-chromosome fit:")
print(pd.DataFrame(best).round(4).to_string(index=False))
"""

md_ub_mix = """## Long-read panel against the mixed panel (SNPs)

Solid: long-read panel (78 founders, no private k-mer filter). Dashed: mixed production panel (231 founders, private
k-mer filter). Same design, units run on both panels. **Not a test of the private k-mer filter**: the two runs also differ
in the founders the pools are drawn from (78 against 231), the number of candidate founders the fit chooses among and
the variant set scored (about 0.86 against 1.03 million Chr1 SNPs). The clean control for the filter is the long-read
panel with the filter on (`benchmarks/p80/data/kmer_pa_p78_filt2inv`), not yet run."""

code_ub_mix = r"""
UM = ["chrom", "bp10k", "bp50k", "ld"]
with PS.style():
    fig, axs = plt.subplots(1, 2, figsize=(PS.FULL_W, 62 * MM), sharey=True)
    for c, (mt, mlab) in enumerate(MAT):
        a = axs[c]
        for panel_, df, ls in (("long-read", UB, "-"), ("mixed", UBM, (0, (4, 2)))):
            for u in UM:
                m = df[(df.cls == "SNP") & (df.mating == mt) & (df.unit == u)].groupby("generations").RMSE.mean()
                a.plot([xg[g] for g in m.index], m.values, color=UCOL[u], ls=ls, lw=PS.LW, marker="o", ms=2)
        a.set_yscale("log"); a.set_xticks(range(len(GENS))); a.set_xticklabels([str(g) for g in GENS]); a.set_xlim(-0.2, len(GENS) - 0.8)
        a.annotate(mlab, (0.03, 0.95), xycoords="axes fraction", ha="left", va="top"); a.set_xlabel("Generations")
    axs[0].set_ylabel("Allele-frequency error\n(RMSE)")
    from matplotlib.lines import Line2D
    hu = [Line2D([], [], color=UCOL[u], lw=PS.LW, label=ULAB[u]) for u in UM]
    hp = [Line2D([], [], color="black", lw=PS.LW, ls=ls, label=f"{p} panel") for p, ls in (("Long-read", "-"), ("Mixed", (0, (4, 2))))]
    fig.legend(handles=hu + hp, loc="lower center", ncol=3, frameon=False, handlelength=2.2, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.12, right=0.98, top=0.97, bottom=0.33, wspace=0.06)
    fig.savefig(f"{PLOTS}/kmate_unitbench_panels.png", dpi=200); plt.show()
t = pd.concat({"long-read": mean_tab(UB, "RMSE")[UM], "mixed": mean_tab(UBM, "RMSE")[UM]}, axis=1)
print("SNP RMSE, mean of 3 seeds, long-read vs mixed panel:"); print(t.round(4).to_string())
"""

nb = nbf.v4.new_notebook()
nb.cells = [nbf.v4.new_markdown_cell(md_title), nbf.v4.new_code_cell(code_load),
            nbf.v4.new_markdown_cell(md_fig), nbf.v4.new_code_cell(code_fig),
            nbf.v4.new_markdown_cell(md_supp),
            nbf.v4.new_markdown_cell(md_s_af), nbf.v4.new_code_cell(code_s_af),
            nbf.v4.new_markdown_cell(md_s_hf), nbf.v4.new_code_cell(code_s_hf),
            nbf.v4.new_markdown_cell(md_s_ft), nbf.v4.new_code_cell(code_s_ft),
            nbf.v4.new_markdown_cell(md_s_lr), nbf.v4.new_code_cell(code_s_lr),
            nbf.v4.new_markdown_cell(md_ub), nbf.v4.new_code_cell(code_ub_load),
            nbf.v4.new_markdown_cell(md_ub_gen), nbf.v4.new_code_cell(code_ub_gen),
            nbf.v4.new_markdown_cell(md_ub_sweep), nbf.v4.new_code_cell(code_ub_sweep),
            nbf.v4.new_markdown_cell(md_ub_mix), nbf.v4.new_code_cell(code_ub_mix)]
nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
nbf.write(nb, NB)
subprocess.run([sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", "--inplace",
                "--ExecutePreprocessor.timeout=600", NB], check=True, cwd=HERE)
print(f"[built+executed] {NB}")
