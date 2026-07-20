#!/usr/bin/env python
"""Analyze beta-binomial+LF GEA: GIF 2x2 table, GIF-calibration, hit-calling,
SNP-vs-nonSNP block comparison, gene annotation, and the figure.

Run in `basic` env (statsmodels/scipy/matplotlib). Writes to betabinom/.
"""
from __future__ import annotations
import sys, os, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
from scipy.stats import chi2, spearmanr
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BB = f"{lib.GEA}/gea_newpanel/results/betabinom"
BIN = f"{lib.GEA}/gea_newpanel/results/binomial_latent"
CM = f"{lib.GEA}/phase1_replication/results/class_matrices"
NULLMED = chi2.ppf(0.5, 1)


def gif(p):
    p = np.asarray(p, float); p = p[np.isfinite(p) & (p > 0)]
    return float(np.median(chi2.isf(p, 1)) / NULLMED)


def gif_calibrate(p, lam):
    p = np.asarray(p, float)
    x2 = chi2.isf(np.clip(p, 1e-300, 1), 1)
    return chi2.sf(x2 / lam, 1)


def load(path):
    d = pd.read_csv(path)
    d = d[np.isfinite(d.pval) & (d.pval > 0)].copy()
    return d


def binom_k0_gif(cls, nsub=20000, seed=3):
    """Raw plain-binomial (K=0, no overdispersion) GIF on a subsample."""
    import statsmodels.api as sm
    from sklearn.preprocessing import StandardScaler
    recs = pd.read_csv(f"{CM}/{cls}_gen9.records.csv")
    pools = pd.read_csv(f"{CM}/gen9.pools.csv")
    clim = pools["bio1"].to_numpy(float)
    N = np.round(pools["total_flowers"].to_numpy(float) * 2).astype(float)
    af = np.load(f"{CM}/{cls}_gen9_af.npy")
    z = StandardScaler().fit_transform(clim.reshape(-1, 1))
    X = np.column_stack([np.ones(len(clim)), z.ravel()])
    rng = np.random.default_rng(seed)
    idx = rng.choice(af.shape[1], min(nsub, af.shape[1]), replace=False)
    pv = []
    for i in idx:
        a = af[:, i]; fin = np.isfinite(a)
        s = np.where(fin, a * N, 0.0); f = np.where(fin, (1 - a) * N, 0.0)
        n = s + f; ok = n > 0
        if ok.sum() < 3 or s[ok].sum() == 0 or f[ok].sum() == 0:
            continue
        try:
            r = sm.GLM(np.column_stack([s[ok], f[ok]]), X[ok],
                       family=sm.families.Binomial()).fit()
            pv.append(r.pvalues[1])
        except Exception:
            pass
    return gif(pv)


def block_peak(d, pcol):
    g = d.groupby("block")[pcol].min().reset_index()
    g["nlp"] = -np.log10(g[pcol].clip(1e-300))
    return g[["block", "nlp"]]


def main():
    genes = lib.load_genes()
    summary = {}
    tables = {}

    for cls in ["snp", "nonsnp"]:
        bb16 = load(f"{BB}/betabinom_lf16_{cls}_gen9_bio1.csv")
        bb0 = load(f"{BB}/betabinom_lf0_{cls}_gen9_bio1.csv")
        b16 = load(f"{BIN}/binomial_lf16_{cls}_gen9_bio1.csv")
        rho = json.load(open(f"{BB}/rho_{cls}.json"))

        g_b0 = binom_k0_gif(cls)
        g_b16 = gif(b16.pval)
        g_bb0 = gif(bb0.pval)
        g_bb16 = gif(bb16.pval)

        # GIF-calibrate the primary model (BB+K16), then call hits on calibrated p
        bb16 = bb16.copy()
        bb16["pcal"] = gif_calibrate(bb16.pval.to_numpy(), g_bb16)
        m = bb16.shape[0]
        bonf = 0.05 / m
        rej, q, _, _ = multipletests(bb16.pcal.to_numpy(), method="fdr_bh")
        bb16["qval"] = q
        n_fdr = int((q < 0.05).sum())
        n_bonf = int((bb16.pcal < bonf).sum())

        summary[cls] = dict(
            rho_global=rho["global"], rho_per_chrom=rho["per_chrom"],
            n_variants=int(m),
            gif_binom_k0=g_b0, gif_binom_k16=g_b16,
            gif_bb_k0=g_bb0, gif_bb_k16=g_bb16,
            n_fdr_q05=n_fdr, n_bonf=n_bonf, bonf_thresh=bonf,
            minp_bb16=float(bb16.pval.min()), minp_cal=float(bb16.pcal.min()),
        )
        tables[cls] = dict(bb16=bb16, b16=b16, bpeak=block_peak(bb16, "pcal"))

        # top-25 by calibrated p, gene annotated
        top = bb16.nsmallest(25, "pcal").copy()
        top = lib.annotate_svs(top, flank=2000, genes=genes)
        cols = ["chrom", "pos", "ref_len", "alt_len", "MAF", "block",
                "slope", "pval", "pcal", "qval", "gene", "gene_name", "n_genes"]
        top[cols].to_csv(f"{BB}/top25_betabinom_lf16_{cls}.csv", index=False)
        print(f"\n=== {cls}: top-8 beta-binom+LF16 hits (GIF-cal) ===")
        print(top[["chrom", "pos", "MAF", "slope", "pcal", "qval",
                   "gene_name", "gene"]].head(8).to_string(index=False))

    # ---- SNP vs non-SNP block-peak comparison ----
    ps = tables["snp"]["bpeak"].rename(columns={"nlp": "nlp_snp"})
    pn = tables["nonsnp"]["bpeak"].rename(columns={"nlp": "nlp_nonsnp"})
    mrg = ps.merge(pn, on="block", how="inner")
    rho_bp, p_bp = spearmanr(mrg.nlp_snp, mrg.nlp_nonsnp)
    summary["block_peak_spearman"] = dict(rho=float(rho_bp), p=float(p_bp),
                                          n_blocks=int(len(mrg)))
    # class-specific top blocks (strong in one, weak in other)
    mrg["diff_nonsnp_minus_snp"] = mrg.nlp_nonsnp - mrg.nlp_snp
    nonsnp_specific = mrg.nlargest(15, "diff_nonsnp_minus_snp")
    snp_specific = mrg.nsmallest(15, "diff_nonsnp_minus_snp")
    nonsnp_specific.to_csv(f"{BB}/nonsnp_specific_blocks.csv", index=False)
    snp_specific.to_csv(f"{BB}/snp_specific_blocks.csv", index=False)

    # ---- Chr1:13,979,249 neighborhood check (non-SNP signal SNPs miss?) ----
    target = 13_979_249
    near = {}
    for cls in ["snp", "nonsnp"]:
        d = tables[cls]["bb16"]
        w = d[(d.chrom == "Chr1") & (d.pos.between(target - 20000, target + 20000))]
        if len(w):
            best = w.nsmallest(1, "pcal").iloc[0]
            near[cls] = dict(n=int(len(w)), best_pos=int(best.pos),
                             best_pcal=float(best.pcal), best_slope=float(best.slope),
                             best_maf=float(best.MAF))
        else:
            near[cls] = dict(n=0)
    summary["chr1_13979249"] = near

    json.dump(summary, open(f"{BB}/betabinom_summary.json", "w"), indent=2, default=str)
    print("\n===== SUMMARY =====")
    print(json.dumps(summary, indent=2, default=str))

    # ---------------------------------------------------------------- FIGURE
    fig = plt.figure(figsize=(14, 10))
    # QQ: plain-binom K16 vs BB K16, both classes
    axq = fig.add_subplot(2, 2, 1)
    for cls, col in [("snp", "#1f77b4"), ("nonsnp", "#d62728")]:
        for tag, dd, ls in [("binom+LF16", tables[cls]["b16"].pval, ":"),
                            ("BB+LF16", tables[cls]["bb16"].pval, "-")]:
            p = np.sort(np.asarray(dd, float))
            p = p[np.isfinite(p) & (p > 0)]
            n = len(p)
            exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
            obs = -np.log10(p)
            # thin for plotting
            k = max(1, n // 6000)
            axq.plot(exp[::k], obs[::k], ls, color=col, lw=1.2,
                     label=f"{cls} {tag} (GIF={gif(dd):.1f})")
    lim = axq.get_xlim()
    axq.plot([0, lim[1]], [0, lim[1]], "k--", lw=0.7)
    axq.set_xlabel("expected -log10 p"); axq.set_ylabel("observed -log10 p")
    axq.set_title("QQ: plain-binomial vs beta-binomial (+K16 LF)")
    axq.legend(fontsize=7)

    # Manhattan, 2 panels (BB K16, GIF-calibrated)
    chrom_order = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
    palette = ["#3b5b92", "#7aa0c4"]
    for pi, cls in enumerate(["snp", "nonsnp"]):
        ax = fig.add_subplot(2, 2, 3 + pi)
        d = tables[cls]["bb16"].copy()
        d = d[d.chrom.isin(chrom_order)]
        d["nlp"] = -np.log10(d.pcal.clip(1e-300))
        off = 0; ticks = []; tlab = []
        for ci, ch in enumerate(chrom_order):
            dc = d[d.chrom == ch].sort_values("pos")
            if not len(dc): continue
            x = off + dc.pos.to_numpy()
            k = max(1, len(dc) // 40000)
            ax.scatter(x[::k], dc.nlp.to_numpy()[::k], s=2,
                       color=palette[ci % 2], rasterized=True)
            ticks.append(off + dc.pos.median()); tlab.append(ch)
            off = x.max() + 1e6
        bonf = summary[cls]["bonf_thresh"]
        ax.axhline(-np.log10(bonf), color="red", lw=0.8, ls="--",
                   label=f"Bonferroni ({summary[cls]['n_bonf']} hits)")
        ax.set_xticks(ticks); ax.set_xticklabels(tlab, fontsize=7)
        ax.set_ylabel("-log10 p (GIF-cal)")
        ax.set_title(f"beta-binomial+K16 GEA — {cls} "
                     f"(GIF {summary[cls]['gif_bb_k16']:.2f}, "
                     f"FDR q<.05: {summary[cls]['n_fdr_q05']})", fontsize=9)
        ax.legend(fontsize=7, loc="upper right")

    fig.suptitle("Overdispersion-aware (beta-binomial) + K=16 latent-factor "
                 "climate-GEA (bio1, gen9)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(f"{BB}/betabinom_qq_manhattan.png", dpi=140)
    print(f"\n-> figure {BB}/betabinom_qq_manhattan.png")
    print(f"-> summary {BB}/betabinom_summary.json")


if __name__ == "__main__":
    main()
