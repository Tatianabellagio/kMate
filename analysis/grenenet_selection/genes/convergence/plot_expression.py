#!/usr/bin/env python
"""Expression in carriers vs non-carriers of a candidate variant, with its lineage control.

Two panels per gene, written as results/plots/expr/expr_<sym>.{png,pdf} and linked into
genes/figures/ as <sym>__expr by organize_figures.py:

  A  1001T rosette expression (raw counts, log scale) for founders that carry the variant
     vs called non-carriers -- the same founders and the same split-record merge as
     carrier_expression_test.py. Each dot is one accession (replicates averaged).
  B  the lineage control: the SAME carrier split tested on 300 random expressed genes.
     Grey = their -log10 p; the accent line is this gene. Carriers can share ancestry,
     and ancestry shifts expression genome-wide, so a difference only counts if it sits
     outside this cloud.

Numbers (n, medians, fold, p, empirical p) are printed to the terminal, not drawn -- the
panels carry no numeric captions by design.

Usage: PY plot_expression.py SYM ...      (symbols from results/carrier_expression_test.csv)
env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd, pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plot_theme as TH                                          # noqa: E402
import carrier_expression_test as CT                             # noqa: E402

TH.apply()
OUT = f"{HERE}/results"
FIG = f"{OUT}/plots/expr"
PROJ = CT.PROJ


def main(syms):
    D = pd.read_csv(f"{OUT}/functional_track_candidates.csv").rename(columns={"symbol": "sym"})
    Gx = pd.read_csv(f"{CT.EQ}/expr_genes.csv")
    Cx = pd.read_csv(f"{CT.EQ}/expr_control.csv")
    graw = Gx.pivot_table(index="acc", columns="id", values="raw", aggfunc="mean")
    ctl = Cx.pivot_table(index="acc", columns="id", values="value", aggfunc="mean")
    os.makedirs(FIG, exist_ok=True)
    for sym in syms:
        v = D[D.sym == sym].iloc[0]
        f = f"{PROJ}/panel/arch3/{v.chrom.lower()}/merged_231_{v.chrom.lower()}_final.vcf.gz"
        car, ref = CT.carriers(pysam.VariantFile(f), v.chrom, int(v.pos), int(v.ref_len),
                               int(v.alt_len))
        r = graw[v.gene]
        cc = [a for a in r.index if a in car and np.isfinite(r[a])]
        rr = [a for a in r.index if a in ref and np.isfinite(r[a])]
        yc, yr = r.reindex(cc).values + 1, r.reindex(rr).values + 1
        p = CT.mw(yc, yr)
        pc = np.array([CT.mw(ctl[c].reindex(cc).values, ctl[c].reindex(rr).values)
                       for c in ctl.columns])
        pc = pc[np.isfinite(pc)]

        fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.2, 4.0),
                                     gridspec_kw=dict(width_ratios=[1.2, 1], wspace=0.35))
        rng = np.random.default_rng(0)
        for x, y, col, lab in [(0, yr, TH.FAINT, "non-carriers"), (1, yc, TH.ACCENT, "carriers")]:
            a1.boxplot([y], positions=[x], widths=0.45, showfliers=False,
                       medianprops=dict(color=TH.LABEL, lw=1.4),
                       boxprops=dict(color=TH.MUTED), whiskerprops=dict(color=TH.MUTED),
                       capprops=dict(color=TH.MUTED))
            a1.scatter(x + rng.uniform(-0.14, 0.14, len(y)), y, s=16, color=col,
                       alpha=0.75, edgecolor="none", zorder=3)
        a1.set_yscale("log")
        # label 1-2-5 steps: GPX6's range spans less than a decade, so decade-only ticks
        # left a single labelled tick on the axis
        from matplotlib.ticker import LogLocator, FuncFormatter
        a1.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
        a1.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        a1.yaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
        a1.set_xticks([0, 1])
        a1.set_xticklabels([f"non-carriers\n(n={len(yr)})", f"carriers\n(n={len(yc)})"])
        a1.set_ylabel(f"{sym} expression (raw counts + 1)")
        a1.set_xlim(-0.6, 1.6)
        TH.panel(a1, "A", f"{v.gene}, 1001T rosettes", y=1.10)

        lp = -np.log10(pc)
        a2.hist(lp, bins=25, color=TH.FAINT, edgecolor="white")
        obs = -np.log10(max(p, 1e-300))
        a2.axvline(obs, color=TH.ACCENT, lw=2)
        a2.set_xlabel("$-\\log_{10}p$, same carrier split")
        a2.set_ylabel("random genes")
        TH.grid_only(a2, "y")
        TH.panel(a2, "B", "lineage control (300 random genes)", y=1.10)
        TH.note(a2, sym, x=min(0.98, obs / max(lp.max(), obs) * 0.95), y=0.95,
                color=TH.ACCENT, ha="right")

        for ext in ("png", "pdf"):
            fig.savefig(f"{FIG}/expr_{sym}.{ext}", dpi=170, bbox_inches="tight")
        plt.close(fig)
        emp = (np.sum(pc <= p) + 1) / (len(pc) + 1)
        print(f"  {sym:<10} carriers {len(yc):>3} median {np.median(yc)-1:7.1f} | non-carriers "
              f"{len(yr):>3} median {np.median(yr)-1:7.1f} | fold {np.median(yc)/np.median(yr):.2f}"
              f" | p = {p:.2g} | empirical p vs lineage control = {emp:.3f}")


if __name__ == "__main__":
    main(sys.argv[1:])
