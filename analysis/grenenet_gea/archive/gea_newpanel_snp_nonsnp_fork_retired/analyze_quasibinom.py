#!/usr/bin/env python
"""Analyze quasibinom outputs: GIF table, QQ + Manhattan figures, gene-annotated
top hits, SNP-vs-non-SNP block-peak comparison. Run in the `basic` env (matplotlib).

Writes to analysis/grenenet_gea/gea_newpanel/results/quasibinom/.
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/gea_newpanel/results/quasibinom"
CHR1_HIT = 13979249      # non-SNP-of-interest locus to check

PCOLS = {"binom": "pval_binom", "quasi": "pval_quasi", "effN": "pval_effN"}


def gif(p):
    p = np.asarray(p, float); p = p[np.isfinite(p) & (p > 0)]
    return np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1)


def load(cls, k=16):
    return pd.read_csv(f"{OUT}/quasibinom_lf{k}_{cls}_gen9_bio1.csv")


def qq_points(p, maxn=40000):
    p = np.sort(np.asarray(p, float)[np.isfinite(p)])
    n = len(p)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    obs = -np.log10(np.clip(p, 1e-320, 1))
    if n > maxn:                       # thin for plotting, keep the tail dense
        idx = np.unique(np.r_[np.linspace(0, n - 1, maxn).astype(int), np.arange(n - 200, n)])
        exp, obs = exp[idx], obs[idx]
    return exp, obs


def main():
    dfs = {cls: load(cls, 16) for cls in ["snp", "nonsnp"]}
    dfs0 = {}
    for cls in ["snp", "nonsnp"]:
        f0 = f"{OUT}/quasibinom_lf0_{cls}_gen9_bio1.csv"
        if os.path.exists(f0):
            dfs0[cls] = pd.read_csv(f0)

    # ---------- 1. GIF table ----------
    rows = []
    for cls in ["snp", "nonsnp"]:
        d = dfs[cls]
        row = {"class": cls, "K": 16, "n_tested": int(np.isfinite(d.pval_binom).sum()),
               "phi_med": float(np.nanmedian(d.phi))}
        for name, col in PCOLS.items():
            row[f"GIF_{name}"] = round(gif(d[col]), 2)
        rows.append(row)
        if cls in dfs0:
            d0 = dfs0[cls]
            row0 = {"class": cls, "K": 0, "n_tested": int(np.isfinite(d0.pval_binom).sum()),
                    "phi_med": float(np.nanmedian(d0.phi))}
            for name, col in PCOLS.items():
                row0[f"GIF_{name}"] = round(gif(d0[col]), 2)
            rows.append(row0)
    gtab = pd.DataFrame(rows)
    gtab.to_csv(f"{OUT}/gif_table.csv", index=False)
    print("=== GIF table ===\n", gtab.to_string(index=False))

    # pick best-calibrated correction (GIF closest to 1 from above, tail survives)
    best = {}
    for cls in ["snp", "nonsnp"]:
        cand = {name: gif(dfs[cls][col]) for name, col in PCOLS.items() if name != "binom"}
        # prefer the one with GIF nearest 1 but >=0.9 (not deflated below)
        best[cls] = min(cand, key=lambda n: abs(cand[n] - 1.0))
    print("best-calibrated per class:", best)

    # ---------- 2. Hit counts (per class, per correction) ----------
    hitrows = []
    for cls in ["snp", "nonsnp"]:
        d = dfs[cls]; ntest = int(np.isfinite(d.pval_binom).sum())
        bonf = 0.05 / ntest
        for name, col in PCOLS.items():
            p = d[col].to_numpy(float); fin = np.isfinite(p)
            q = np.full(len(p), np.nan); q[fin] = lib.bh(p[fin])
            hitrows.append({"class": cls, "correction": name,
                            "n_bonf": int((p[fin] < bonf).sum()),
                            "n_fdr05": int((q[fin] < 0.05).sum()),
                            "bonf_thresh": bonf})
    htab = pd.DataFrame(hitrows)
    htab.to_csv(f"{OUT}/hit_counts.csv", index=False)
    print("\n=== hit counts ===\n", htab.to_string(index=False))

    # ---------- 3. QQ figure ----------
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    colr = {"binom": "#999999", "quasi": "#1f77b4", "effN": "#d62728"}
    lbl = {"binom": "binomial+LF16", "quasi": "quasi-binom+LF16", "effN": "effective-N+LF16"}
    for ax, cls in zip(axes, ["snp", "nonsnp"]):
        d = dfs[cls]
        mx = 0
        for name, col in PCOLS.items():
            ex, ob = qq_points(d[col])
            ax.scatter(ex, ob, s=4, c=colr[name], alpha=0.5,
                       label=f"{lbl[name]} (GIF {gif(d[col]):.2f})")
            mx = max(mx, ex.max(), ob.max())
        ax.plot([0, mx], [0, mx], "k--", lw=0.8)
        ax.set_title(f"{cls}  (n={int(np.isfinite(d.pval_binom).sum()):,})")
        ax.set_xlabel("expected -log10(p)"); ax.set_ylabel("observed -log10(p)")
        ax.legend(fontsize=7, loc="upper left")
    fig.suptitle("Overdispersion correction QQ: binomial vs quasi-binomial vs effective-N (all +K16)")
    fig.tight_layout()
    fig.savefig(f"{OUT}/qq_quasibinom.png", dpi=140)
    plt.close(fig)
    print(f"\nwrote {OUT}/qq_quasibinom.png")

    # ---------- 4. Manhattan (best-calibrated) ----------
    chrom_order = [f"Chr{i}" for i in range(1, 6)]
    fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=False)
    for ax, cls in zip(axes, ["snp", "nonsnp"]):
        d = dfs[cls].copy()
        col = PCOLS[best[cls]]
        d = d[np.isfinite(d[col]) & (d.chrom.isin(chrom_order))].copy()
        d["mlp"] = -np.log10(np.clip(d[col].to_numpy(float), 1e-320, 1))
        # cumulative x
        off = 0; ticks = []; xall = np.zeros(len(d));
        cvals = d.chrom.to_numpy(); pvals = d.pos.to_numpy()
        for ci, c in enumerate(chrom_order):
            m = cvals == c
            xall[m] = pvals[m] + off
            if m.any():
                ticks.append((off + pvals[m].max() / 2, c))
                off += pvals[m].max() + 5e6
        col_cyc = np.where(np.isin(cvals, chrom_order[::2]), "#2c7fb8", "#7fcdbb")
        ax.scatter(xall, d.mlp, s=3, c=col_cyc, alpha=0.6)
        ntest = int(np.isfinite(dfs[cls].pval_binom).sum())
        ax.axhline(-np.log10(0.05 / ntest), color="red", lw=0.8, ls="--",
                   label=f"Bonferroni 0.05/{ntest:,}")
        # FDR line: p at largest q<0.05
        p = dfs[cls][col].to_numpy(float); fin = np.isfinite(p)
        q = lib.bh(p[fin]); pf = p[fin]
        if (q < 0.05).any():
            fdr_p = pf[q < 0.05].max()
            ax.axhline(-np.log10(fdr_p), color="orange", lw=0.8, ls=":",
                       label=f"FDR q<0.05 (p<{fdr_p:.1e})")
        ax.set_xticks([t[0] for t in ticks]); ax.set_xticklabels([t[1] for t in ticks])
        ax.set_ylabel("-log10(p)")
        ax.set_title(f"{cls}  {best[cls]}+LF16  (GIF {gif(dfs[cls][col]):.2f})")
        ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(f"{OUT}/manhattan_best.png", dpi=140)
    plt.close(fig)
    print(f"wrote {OUT}/manhattan_best.png")

    # ---------- 5. gene-annotated top-25 per class ----------
    genes = lib.load_genes()
    for cls in ["snp", "nonsnp"]:
        d = dfs[cls].copy()
        col = PCOLS[best[cls]]
        d = d[np.isfinite(d[col])].copy()
        top = d.nsmallest(25, col)
        top = lib.annotate_svs(top, flank=2000, genes=genes)
        keep = ["chrom", "pos", "ref_len", "alt_len", "MAF", "block", "slope",
                "phi", "pval_binom", "pval_quasi", "pval_effN", "gene", "gene_name", "genes_all"]
        top[keep].to_csv(f"{OUT}/top25_{cls}.csv", index=False)
        print(f"\n=== top-25 {cls} by {best[cls]} ===")
        print(top[["chrom", "pos", "MAF", "slope", "phi", col, "gene_name", "gene"]].to_string(index=False))

    # ---------- 6. SNP vs non-SNP block-peak comparison ----------
    peaks = {}
    for cls in ["snp", "nonsnp"]:
        d = dfs[cls]; col = PCOLS[best[cls]]
        dd = d[np.isfinite(d[col]) & (d.block.astype(str) != "")].copy()
        pk = dd.groupby("block")[col].min().rename(f"p_{cls}")
        peaks[cls] = pk
    merged = pd.concat([peaks["snp"], peaks["nonsnp"]], axis=1, join="inner").dropna()
    rho, prho = stats.spearmanr(-np.log10(merged.p_snp.clip(1e-320)),
                                -np.log10(merged.p_nonsnp.clip(1e-320)))
    print(f"\n=== block-peak SNP vs non-SNP ===\n"
          f"common blocks={len(merged)}  Spearman rho(-log10 p_peak)={rho:.3f} (p={prho:.1e})")

    # class-specific top blocks (strong in one, weak in other)
    m = merged.copy()
    m["mlp_snp"] = -np.log10(m.p_snp.clip(1e-320)); m["mlp_nonsnp"] = -np.log10(m.p_nonsnp.clip(1e-320))
    m["diff"] = m.mlp_nonsnp - m.mlp_snp
    nonsnp_specific = m.nlargest(15, "diff")
    snp_specific = m.nsmallest(15, "diff")
    m.to_csv(f"{OUT}/block_peaks_snp_vs_nonsnp.csv")
    nonsnp_specific.to_csv(f"{OUT}/nonsnp_specific_blocks.csv")
    snp_specific.to_csv(f"{OUT}/snp_specific_blocks.csv")

    # annotate class-specific top blocks with genes (use each block's lead variant)
    def annotate_blocks(blocks_idx, cls):
        d = dfs[cls]; col = PCOLS[best[cls]]
        sub = d[d.block.isin(blocks_idx)].copy()
        lead = sub.loc[sub.groupby("block")[col].idxmin()]
        lead = lib.annotate_svs(lead, flank=2000, genes=genes)
        return lead[["block", "chrom", "pos", "MAF", "slope", col, "gene_name", "gene"]]
    print("\n--- non-SNP-specific top blocks (genes, via non-SNP lead) ---")
    print(annotate_blocks(nonsnp_specific.index, "nonsnp").to_string(index=False))
    print("\n--- SNP-specific top blocks (genes, via SNP lead) ---")
    print(annotate_blocks(snp_specific.index, "snp").to_string(index=False))

    # ---------- 7. Chr1 locus of interest ----------
    print(f"\n=== near Chr1:{CHR1_HIT:,} (+/-5kb) ===")
    for cls in ["snp", "nonsnp"]:
        d = dfs[cls]
        near = d[(d.chrom == "Chr1") & (d.pos.between(CHR1_HIT - 5000, CHR1_HIT + 5000))]
        if len(near):
            col = PCOLS[best[cls]]
            best_near = near.loc[near[col].idxmin()]
            print(f"  {cls}: {len(near)} variants; best pos={int(best_near.pos)} "
                  f"MAF={best_near.MAF:.3f} slope={best_near.slope:.3f} "
                  f"p_binom={best_near.pval_binom:.1e} p_quasi={best_near.pval_quasi:.1e} "
                  f"p_effN={best_near.pval_effN:.1e}")
        else:
            print(f"  {cls}: no variants in window")

    print("\nDONE")


if __name__ == "__main__":
    main()
