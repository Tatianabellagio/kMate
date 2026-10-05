#!/usr/bin/env python
"""Build + execute notebooks/candidates_v2_<NAME>.ipynb -- locus figures for genes of the NEW candidate list.

Load-only: reads results/candidate_list_v2_events_cr50/locus/<GENE>/ (written by v2/locus_data.py) and v2_candidate_list.csv.
Per gene: (1) locus view per significant climate axis -- LFMM association of non-SNP units (green) and SNPs (grey) along the chromosome,
Bonferroni lines, the candidate unit(s) marked, gene models and the ATAC / TE tracks; (2) per-pool frequency change (delta p, gen 9 - founders)
against the climate axis, with the OLS line; (3) carrier vs non-carrier expression of the target gene. No chart titles (project convention).
usage: python _build_locus_notebook.py NAME GENE [GENE ...]      env: basic
"""
import os, sys
import pandas as pd
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
GS = os.path.abspath(os.path.join(CONV, "..", ".."))
V2 = f"{CONV}/results/candidate_list_v2_events_cr50"
NAME, GENES = sys.argv[1], sys.argv[2:]
OUT = f"{GS}/notebooks/candidates_v2_{NAME}.ipynb"

L = pd.read_csv(f"{V2}/v2_candidate_list.csv").set_index("gene")

setup = r"""
import os, sys, glob
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
GS = "__GS__"; V2 = "__V2__"; PLOTS = f"{GS}/notebooks/plots"
sys.path.insert(0, f"{GS}/genes/convergence")
import plot_theme as TH
TH.apply(); plt.rcParams["figure.dpi"] = 120
SYM = pd.read_csv(f"{V2}/v2_candidate_list.csv").set_index("gene").symbol.to_dict()

def load(g):
    d = f"{V2}/locus/{g}"
    R = dict(units=pd.read_csv(f"{d}/units.csv", low_memory=False), genes=pd.read_csv(f"{d}/genes_window.csv"),
             te=pd.read_csv(f"{d}/te_window.csv"), atac=pd.read_csv(f"{d}/atac_window.csv"), dp=pd.read_csv(f"{d}/dp_pools.csv"),
             zoom={os.path.basename(f)[5:-4]: pd.read_csv(f) for f in sorted(glob.glob(f"{d}/zoom_*.csv"))})
    R["expr"] = {int(os.path.basename(f)[5:-4]): pd.read_csv(f) for f in sorted(glob.glob(f"{d}/expr_*.csv"))}
    return R

def pack_genes(G, lo, hi, g):
    # greedy row packing so gene names never collide: interval = gene body widened to its label
    genes = sorted(G[G.feature == "gene"].itertuples(), key=lambda r: r.start); span = hi - lo; rows, placed = [], []
    for gr in genes:
        name = SYM.get(gr.gene, gr.gene); name = name if isinstance(name, str) else gr.gene
        cx = (gr.start + gr.end) / 2; w = span * 0.0068 * len(name)
        a, b = min(gr.start, cx - w / 2), max(gr.end, cx + w / 2)
        lvl = 0
        while any(not (b + span * 0.004 < pa or a - span * 0.004 > pb) for (pl, pa, pb) in placed if pl == lvl):
            lvl += 1
        placed.append((lvl, a, b)); rows.append((gr, name, lvl))
    return rows

def locus_fig(g, axis, R, slug):
    U = R["units"]; Z = R["zoom"][axis]; ch = U.chrom.iat[0]
    lo, hi = Z.pos.min(), Z.pos.max(); G = R["genes"]; rows = pack_genes(G, lo, hi, g); nlev = max([l for _, _, l in rows] + [0]) + 1
    fig, ax = plt.subplots(4, 1, figsize=(11, 6.4 + 0.45 * nlev), sharex=True, gridspec_kw=dict(height_ratios=[4, 0.7 + 0.55 * nlev, 0.4, 0.4], hspace=0.08))
    a = ax[0]; yl = []
    for cls, col, sz, lab in (("snp", TH.FAINT, 3, "SNPs"), ("nonsnp", TH.GREENS[1], 9, "non-SNP units (indels, SV events)")):
        d = Z[Z.cls == cls]; a.scatter(d.pos / 1e3, -np.log10(d.pval.clip(lower=1e-300)), s=sz, color=col, alpha=0.7 if cls == "snp" else 0.85, linewidths=0, label=lab, zorder=2)
        yb = -np.log10(0.05 / d.n_tested.iat[0]); yl.append(yb)
        a.axhline(yb, color=col, lw=0.9, ls=(0, (4, 3)), zorder=1, label=f"Bonferroni ({'SNPs' if cls == 'snp' else 'non-SNP units'})")
    pts = []
    for u in U.itertuples():
        row = Z[(Z.cls == "nonsnp") & (Z.unit == u.unit)]
        if len(row):
            pts.append((u.pos, -np.log10(max(row.pval.iat[0], 1e-300)), f"{'+' if u.alt_len > u.ref_len else '−'}{abs(int(u.alt_len - u.ref_len))} bp"))
    for p, y, _ in pts: a.scatter([p / 1e3], [y], s=70, color=TH.ACCENT, edgecolor="white", linewidths=0.8, zorder=4)
    grp = []                                                      # label variants closer than 1.2% of the window together
    for p, y, t in sorted(pts):
        if grp and p - grp[-1][-1][0] < 0.012 * (hi - lo): grp[-1].append((p, y, t))
        else: grp.append([(p, y, t)])
    for gp in grp:
        a.annotate(" / ".join(t for _, _, t in gp), (np.mean([p for p, _, _ in gp]) / 1e3, max(y for _, y, _ in gp)), xytext=(8, 8), textcoords="offset points", fontsize=9, color=TH.ACCENT)
    a.set_ylim(-0.4, max([y for _, y, _ in pts] + yl) * 1.12)
    a.set_ylabel(f"−log10 p ({axis})")
    a.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=4, markerscale=1.8, fontsize=8.5, handlelength=1.8)
    b = ax[1]
    for gr, name, lvl in rows:
        y = -lvl * 1.0; isT = gr.gene == g; col = TH.ACCENT if isT else TH.MUTED
        b.plot([gr.start / 1e3, gr.end / 1e3], [y, y], color=col, lw=1.1, zorder=2)
        for f, h in (("exon", 0.18), ("CDS", 0.34), ("five_prime_UTR", 0.12), ("three_prime_UTR", 0.12)):
            for e in G[(G.gene == gr.gene) & (G.feature == f)].itertuples():
                b.add_patch(Rectangle((e.start / 1e3, y - h / 2), (e.end - e.start) / 1e3, h, color=col, lw=0, zorder=3 if f == "CDS" else 2, alpha=0.95 if f == "CDS" else 0.6))
        tip = gr.end / 1e3 if gr.strand == "+" else gr.start / 1e3; sg = 1 if gr.strand == "+" else -1
        b.annotate("", (tip + sg * 0.6, y), xytext=(tip - sg * 0.6, y), arrowprops=dict(arrowstyle="-|>", color=col, lw=0.8, mutation_scale=7))
        b.annotate(name, ((gr.start + gr.end) / 2e3, y + 0.30), ha="center", fontsize=8, color=TH.ACCENT if isT else TH.TICK)
    b.set_ylim(-nlev + 0.45, 0.85); b.set_yticks([]); b.grid(False); b.set_ylabel("genes")
    c = ax[2]
    for e in R["atac"].itertuples():
        c.add_patch(Rectangle((e.start / 1e3, 0), (e.end - e.start) / 1e3, 1, color=TH.GREENS[2], alpha=min(1.0, 0.25 + 0.17 * e.n_tissues), lw=0))
    c.set_ylim(0, 1); c.set_yticks([]); c.grid(False); c.set_ylabel("ATAC", rotation=0, ha="right", va="center")
    d_ = ax[3]
    for e in R["te"].itertuples():
        d_.add_patch(Rectangle((e.start / 1e3, 0), (e.end - e.start) / 1e3, 1, color=TH.MUTED, alpha=0.7, lw=0))
    d_.set_ylim(0, 1); d_.set_yticks([]); d_.grid(False); d_.set_ylabel("TEs", rotation=0, ha="right", va="center")
    for p, _, _ in pts:
        for q in ax: q.axvline(p / 1e3, color=TH.ACCENT, lw=0.7, ls=(0, (2, 3)), alpha=0.7, zorder=0)
    ax[3].set_xlabel(f"{ch} position (kb)"); ax[3].set_xlim(lo / 1e3, hi / 1e3)
    fig.savefig(f"{PLOTS}/{slug}_locus_{axis}.png", dpi=200, bbox_inches="tight"); plt.show()

def freq_fig(g, R, slug):
    U = R["units"]; D = R["dp"]; pairs = []
    for u in U.itertuples():
        axes = [a for a in str(u.sig_axes).split(",") if a in D.columns] or [u.best_axis]
        pairs += [(u, a) for a in axes]
    fig, axs = plt.subplots(1, len(pairs), figsize=(4.7 * len(pairs), 4.4), squeeze=False)
    for ax, (u, a) in zip(axs[0], pairs):
        x, y = D[a].to_numpy(float), D[f"dp_unit{u.unit}"].to_numpy(float)
        ss = pd.DataFrame({"x": x, "y": y, "s": D.site.to_numpy()}).groupby("s").mean()
        h1 = ax.scatter(x, y, s=10, color=TH.FAINT, linewidths=0, zorder=2, label=f"pools (n={len(D)})")
        h2 = ax.scatter(ss.x, ss.y, s=34, color=TH.GREENS[2], edgecolor="white", linewidths=0.6, zorder=3, label=f"site means (n={len(ss)})")
        k, c0 = np.polyfit(x, y, 1); xs = np.linspace(x.min(), x.max(), 50)
        h3, = ax.plot(xs, k * xs + c0, color=TH.ACCENT, lw=1.4, zorder=4, label="least-squares line over pools"); r = np.corrcoef(x, y)[0, 1]
        ax.axhline(0, color=TH.AXIS, lw=0.7, zorder=1)
        ax.set_xlabel(f"{a} (climate of the garden)"); ax.set_ylabel("allele frequency change Δp (gen 9 − founders)")
        ax.set_ylim(top=ax.get_ylim()[1] + 0.12 * (ax.get_ylim()[1] - ax.get_ylim()[0]))
        TH.note(ax, f"unit {u.unit}: {'+' if u.alt_len > u.ref_len else '−'}{abs(int(u.alt_len - u.ref_len))} bp, MAF {u.MAF:.2f}\nr = {r:.2f} over pools; LFMM −log10 p = {u.gea_nlp:.1f}", y=0.98, x=0.02)
    fig.legend([h1, h2, h3], [h1.get_label(), h2.get_label(), h3.get_label()], loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig(f"{PLOTS}/{slug}_freq_climate.png", dpi=200, bbox_inches="tight"); plt.show()

def expr_fig(g, R, slug):
    U = R["units"].set_index("unit"); n = len(R["expr"])
    fig, axs = plt.subplots(1, n, figsize=(4.0 * n, 4.2), squeeze=False)
    rng = np.random.default_rng(0)
    for ax, (u, E) in zip(axs[0], R["expr"].items()):
        r = U.loc[u]; allv = []
        for i, (k, col) in enumerate((("noncarrier", TH.FAINT), ("carrier", TH.ACCENT))):
            v = np.log2(E[(E.status == k)].raw.dropna().to_numpy() + 1); allv += list(v)
            ax.boxplot([v], positions=[i], widths=0.5, showfliers=False, medianprops=dict(color=TH.LABEL), boxprops=dict(color=TH.MUTED), whiskerprops=dict(color=TH.MUTED), capprops=dict(color=TH.MUTED))
            ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=9, color=col, alpha=0.75, linewidths=0, zorder=3)
        lo_, hi_ = min(allv), max(allv); ax.set_ylim(lo_ - 0.18 * (hi_ - lo_), hi_ + 0.62 * (hi_ - lo_))
        for i, k in enumerate(("noncarrier", "carrier")):
            ax.annotate(f"n = {int((E.status == k).sum() - E[E.status == k].raw.isna().sum())}", (i, ax.get_ylim()[0]), xytext=(0, 5), textcoords="offset points", ha="center", fontsize=8, color=TH.MUTED)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["non-carriers", "carriers"]); ax.set_ylabel(f"{SYM.get(g, g)} expression, log2(count + 1)"); ax.grid(axis="x", visible=False)
        st = r.get("status", "")
        txt = (f"unit {u}: carrier / non-carrier = {r.fold:.2f}\nMann-Whitney p = {r.p_gene:.2g}\nlineage-control p = {r.p_emp:.3f}" if st == "tested" else f"unit {u}: {st}")
        TH.note(ax, txt, y=0.98, x=0.03)
    fig.tight_layout(); fig.savefig(f"{PLOTS}/{slug}_expression.png", dpi=200, bbox_inches="tight"); plt.show()
"""
setup = setup.replace("__GS__", GS).replace("__V2__", V2)

cells = [new_markdown_cell(f"""# Locus figures — {', '.join(GENES)}

Genes of the **new candidate list** (`candidate_list_v2_events_cr50`: GEA only, SVs as Truvari events, MAF ≥ 0.05, call rate ≥ 50%, raw LFMM p; no GWAS).
Each gene: **(1)** the locus on every significant climate axis (non-SNP units green, SNPs grey, dashed lines = per-class Bonferroni, the candidate unit in red; gene models, ATAC peaks, reference TEs);
**(2)** the per-pool allele-frequency change against the garden's climate; **(3)** carrier versus non-carrier expression of the gene (1001 Transcriptomes).
Data: `results/candidate_list_v2_events_cr50/locus/<gene>/` (written by `v2/locus_data.py`). Raw uncalibrated p: a candidate net, not a calibrated hit list."""),
         new_code_cell(setup.strip())]
for g in GENES:
    r = L.loc[g]; sym = r.symbol
    cells.append(new_markdown_cell(
        f"## {sym} ({g}) — {r.themes if isinstance(r.themes, str) else 'no theme'}\n\n"
        f"Rank {int(r['rank'])} of {len(L):,}. Best functional handle: **{r.best_ftier}** ({r.best_vclass}, {int(r.best_size)} bp); "
        f"extra lines of evidence in its best unit: **{int(r.best_n_extra_lines)}** (gated GEA line: {bool(r.L_gea)}); GEA recurrence: {int(r.gea_n_clusters)} climate cluster(s), cleanest axis λ = {r.gea_min_lam:.2f}."))
    cells.append(new_code_cell(f'g = "{g}"; slug = "{sym}"; R = load(g)\nR["units"][["unit","chrom","pos","ref_len","alt_len","vclass","ftier","MAF","gea_best_axis","gea_nlp","gea_n_clusters","in_atac","tfbs_lost","tfbs_gained","best_r2_shortread","fold","p_emp","U_gea","U_expr","U_regulatory","U_snp_blind"]].round(3)'))
    cells.append(new_markdown_cell(f"### {sym}: locus on the significant climate axes"))
    cells.append(new_code_cell('for axis in R["zoom"]:\n    locus_fig(g, axis, R, slug)'))
    cells.append(new_markdown_cell(f"### {sym}: allele-frequency change against climate"))
    cells.append(new_code_cell("freq_fig(g, R, slug)"))
    cells.append(new_markdown_cell(f"### {sym}: carrier vs non-carrier expression"))
    cells.append(new_code_cell("expr_fig(g, R, slug)"))

nb = new_notebook(cells=cells)
ExecutePreprocessor(timeout=900, kernel_name="python3").preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
nbf.write(nb, OUT)
print("[saved]", OUT)
