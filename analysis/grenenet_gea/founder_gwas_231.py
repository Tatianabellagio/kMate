#!/usr/bin/env python
"""Founder-panel fitness GWAS on ALL 231 founders, LINEAR scale + LOCO (default site 4).

Individuals = the 231 founders (all are in the 8-rep seed mix; kMate ~0 estimates are an
ESTIMATION artifact, not real absence -> full panel). Trait = each founder's selection response
= LINEAR-scale slope of its genome-wide frequency over gens 0-3 at the site (matches the linear
block LD-LMM fix: raw-frequency slope, NO logit floor that inflated rare founders), flower-weighted
across replicate plots, then quantile-normalized (QN is load-bearing against the 9764-type winner
dominating the marginal fit). Genotype = founder's haploblock allele (k-1 one-vs-rest per dynld-K500
block). Kinship random effect controls relatedness (EMMAX/P3D).

MAF / GRM decoupling (KEY):
  * TEST down to MAF>=1% (MAC>=MAC_MIN=3).
  * Build the kinship GRM from COMMON markers only (MAF>=5%, MAC>=MAC_GRM=12), LEAVE-ONE-CHROMOSOME-
    OUT. Admitting rare markers into the GRM pushes the variance component delta->0 (h2->1) on some
    chromosomes and AMPLIFIES markers with no marginal signal (p_naive~0.9 -> p_kin~3e-6, pure
    delta-degeneracy artifact). Decoupling GRM(common) from test(all) keeps delta stable.
The MAC distribution of any hits is printed as a standing diagnostic (rare-marker inflation guard).

REJECTED APPROACH -- DO NOT RE-ADD (2026-06-26): a replicate-PRECISION arm (heteroscedastic residual
e_f~N(0, se2 + D_f) weighting founders by inverse replicate variance, even with D_f FLOORED at the
flower-census binomial sampling variance). It does not help and is artifact-prone: (a) the floored
D_f is near-uniform across founders (~87% hit the same floor) so it barely re-weights; (b) on the raw
(non-QN) scale it leaks rare-founder leverage -> manufactures rare-marker false positives (median MAC
~4) at the TEST level regardless of the GRM fix; (c) it is only clean where it finds nothing (0 hits
at MAF>=5%). The naive 1/SE^2 version was worse (h2->0.99, degenerate V, 13 spurious Bonferroni).
Net: the honest model is unweighted + QN; replicate precision adds nothing. See memory
founder-gwas-audit-and-fix.

Writes site<ID>_founder_gwas231.csv + _meta.json + manhattan. Env: kmate. SITE via env.
"""
import os, sys, glob, json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from founder_genotype import build_genotype, emma_reml_delta
from ecotype_selection_site import genome_h

H = "results/grenenet_gea/hapfreq"
WIN = "results/grenenet_kmate_window"
SEED = "results/grenenet_kmate_window_seedmix"
SITE = int(os.environ.get("SITE", 4))
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
Tg = np.array([0.0, 1.0, 2.0, 3.0]); coef = (Tg - 1.5)
MAC_MIN = int(os.environ.get("MAC_MIN", 3))                # TEST set: ~MAF>=1% of 231 founders (MAC 2 = 0.87%)
MAC_GRM = int(os.environ.get("MAC_GRM", 12))               # GRM set: ~MAF>=5% (common only -> stable delta)


def build_trait(nF):
    """Per-founder LINEAR genome-wide selection slope, flower-weighted across replicate plots."""
    seeds = sorted({p.split("/")[-1].split("_Chr")[0]
                    for p in glob.glob(f"{SEED}/*_Chr1.h_blocks_per_chrom.npz")})
    p0 = np.mean([genome_h(s, SEED) for s in seeds], 0)               # founding freqs (231)
    pt = lib.pool_table(); s = pt[pt["site"] == SITE].copy()
    s = s[s.sampleid.astype(str).apply(lambda x: os.path.exists(f"{WIN}/{x}_Chr1.h_blocks_per_chrom.npz"))]
    cellH = {}
    for (gen, plot), g in s.groupby(["generation", "plot"]):
        hs, ws = [], []
        for _, r in g.iterrows():
            h = genome_h(str(r.sampleid), WIN)
            if h is not None:
                w = r.flowerscollected if np.isfinite(r.flowerscollected) and r.flowerscollected > 0 else 1.0
                hs.append(h); ws.append(w)
        if hs:                                                        # flower-weighted cell mean
            ws = np.asarray(ws); cellH[(int(gen), int(plot))] = (np.vstack(hs) * ws[:, None]).sum(0) / ws.sum()
    present = {}
    for (gen, plot) in cellH:
        present.setdefault(plot, set()).add(gen)
    plots = sorted(pl for pl, gs in present.items() if {1, 2, 3} <= gs)
    n = len(plots)
    slopes = np.zeros((n, nF))
    for j, pl in enumerate(plots):
        y = np.vstack([p0, cellH[(1, pl)], cellH[(2, pl)], cellH[(3, pl)]])
        slopes[j] = (coef[:, None] * np.clip(y, 0.0, 1.0)).sum(0) / 5.0
    return slopes.mean(0), n


def main():
    G, founders, reg = build_genotype(); nF = len(founders)
    reg["unit"] = reg.chrom + ":" + reg.start.astype(str) + "-" + reg.end.astype(str)
    s_f, n = build_trait(nF)
    print(f"site {SITE}: {nF} founders, trait = LINEAR genome-wide selection slope over {n} plots")

    # ---- markers: drop per-block reference (k-1); TEST set MAF>=1%, GRM set common MAF>=5% ----
    cnt = G.sum(0); blk = reg.block.to_numpy()
    order = np.lexsort((-cnt, blk)); sbk = blk[order]
    is_ref = np.zeros(len(cnt), bool)
    is_ref[order[np.concatenate([[True], sbk[1:] != sbk[:-1]])]] = True
    n_premaf = int(((~is_ref) & (cnt >= 1) & (cnt <= nF - 1)).sum())    # all polymorphic (pre-MAF)
    poly = (~is_ref) & (cnt >= MAC_MIN) & (cnt <= nF - MAC_MIN)         # TEST set (MAF>=1%)
    grm_set = (~is_ref) & (cnt >= MAC_GRM) & (cnt <= nF - MAC_GRM)      # GRM set (common, MAF>=5%)
    Gp = G[:, poly].astype(np.float64); regp = reg[poly].reset_index(drop=True); M = Gp.shape[1]
    chrom_m = regp.chrom.to_numpy()
    Gg = G[:, grm_set].astype(np.float64); chrom_g = reg.chrom.to_numpy()[grm_set]; Mg = Gg.shape[1]
    print(f"  TEST markers: {M:,} (MAC>={MAC_MIN}, MAF>=~{MAC_MIN/nF:.1%}); GRM from {Mg:,} COMMON "
          f"markers (MAC>={MAC_GRM}, MAF>=~{MAC_GRM/nF:.0%}); {n_premaf:,} polymorphic total")

    # ---- LOCO GRMs from COMMON markers OFF the focal chrom -> stable variance components ----
    def grm(cols):
        sub = Gg[:, cols]; p = sub.mean(0)
        Z = (sub - p) / np.sqrt(p * (1 - p) + 1e-9)
        K = (Z @ Z.T) / cols.sum(); return (K + K.T) / 2

    # ---- trait: quantile-normalized; naive (no-kinship) scan ----
    yq = stats.norm.ppf((stats.rankdata(s_f) - 0.5) / nF)
    yc = yq - yq.mean(); Gc = Gp - Gp.mean(0)
    rr = (yc @ Gc) / np.sqrt((Gc ** 2).sum(0) * (yc @ yc) + 1e-30)
    t_naive = rr * np.sqrt((nF - 2) / np.clip(1 - rr ** 2, 1e-12, None))
    p_naive = 2 * stats.t.sf(np.abs(t_naive), df=nF - 2)

    # ---- kinship-corrected EMMAX (LOCO, P3D per marker) ----
    beta_k = np.full(M, np.nan); p_k = np.full(M, np.nan); delta_by_chr = {}
    for ch in CHROMS:
        on = (chrom_m == ch)
        if not on.any():
            continue
        K = grm(chrom_g != ch)                                       # LOCO GRM from common markers off-chrom
        lam, U = np.linalg.eigh(K); lam = np.clip(lam, 1e-9, None)
        yr = U.T @ yq; Xr = U.T @ np.ones((nF, 1))
        delta = emma_reml_delta(yr, Xr, lam); delta_by_chr[ch] = float(delta)
        w = 1.0 / (lam + delta); Br = U.T @ Gp[:, on]; a = Xr[:, 0]
        Saa = (w * a * a).sum(); Say = (w * a * yr).sum()
        Sab = ((w * a)[None, :] @ Br).ravel(); Sbb = (w[:, None] * Br ** 2).sum(0)
        Sby = ((w * yr)[None, :] @ Br).ravel()
        det = Saa * Sbb - Sab ** 2
        b = (Saa * Sby - Sab * Say) / np.clip(det, 1e-30, None)
        alpha = (Sbb * Say - Sab * Sby) / np.clip(det, 1e-30, None)
        rss = (w * yr * yr).sum() - alpha * Say - b * Sby
        s2 = np.clip(rss, 1e-30, None) / (nF - 2)
        se = np.sqrt(np.clip(s2 * Saa / np.clip(det, 1e-30, None), 1e-30, None))
        beta_k[on] = b; p_k[on] = 2 * stats.t.sf(np.abs(b / se), df=nF - 2)

    def lamgc(pv): return np.median(stats.chi2.isf(np.clip(pv, 1e-300, 1), 1)) / stats.chi2.ppf(0.5, 1)
    def bh(pv):
        m = len(pv); o = pv.argsort(); q = np.empty(m)
        q[o] = np.minimum.accumulate((pv[o] * m / (np.arange(m) + 1))[::-1])[::-1]; return np.clip(q, 0, 1)

    d = regp[["chrom", "start", "end", "unit"]].copy()
    d["mac"] = cnt[poly]; d["beta"] = beta_k; d["p_naive"] = p_naive; d["p_kin"] = p_k; d["q_kin"] = bh(p_k)
    bonf = 0.05 / M
    d.sort_values("p_kin").to_csv(f"{H}/site{SITE}_founder_gwas231.csv", index=False)

    nb = int((d.p_kin < bonf).sum()); nq = int((d.q_kin < 0.05).sum())
    hits = d[d.q_kin < 0.05]
    macs = (f"median MAC {int(hits.mac.median())}, frac MAC<=20 {(hits.mac<=20).mean():.2f}" if nq else "-")
    lgc_n, lgc_k = lamgc(p_naive), lamgc(p_k)
    print(f"\n  naive lambda_GC = {lgc_n:.2f}  ->  kinship-corrected (LOCO) lambda_GC = {lgc_k:.2f}")
    print(f"  LOCO delta by chrom: {({k: round(v,3) for k,v in delta_by_chr.items()})}")
    print(f"  hits: Bonferroni {nb} | FDR q<0.05 {nq} | best q={d.q_kin.min():.3g}")
    print(f"  hit MAC profile: {macs}   (overall median MAC {int(d.mac.median())})")
    print("\n  TOP kinship-corrected blocks (suggestive; none clear FDR if nq=0):")
    print(d.sort_values("p_kin").head(10)[["chrom", "start", "end", "mac", "beta", "p_naive", "p_kin", "q_kin"]].to_string(index=False))

    meta = dict(site=SITE, n_founders=nF, M=int(M), n_plot=n, mac_min=MAC_MIN, maf_min=MAC_MIN / nF,
                mac_grm=MAC_GRM, M_grm=int(Mg), n_dropped_rare=int(n_premaf - M),
                lambda_naive=float(lgc_n), lambda_kin=float(lgc_k), delta_by_chrom=delta_by_chr,
                n_bonf=nb, n_fdr=nq, best_q=float(d.q_kin.min()),
                model="linear-scale founder GWAS, QN trait, LOCO GRM from common markers (MAF>=5%), test to MAF>=1%")
    json.dump(meta, open(f"{H}/site{SITE}_founder_gwas231_meta.json", "w"), indent=2)

    # ---- Manhattan: naive vs kinship-corrected ----
    g2 = d.copy(); g2["mid"] = (g2.start + g2.end) / 2
    g2 = g2.sort_values(["chrom", "start"]).reset_index(drop=True)
    off, centers, x = 0.0, [], np.zeros(len(g2))
    for ch in CHROMS:
        mk = (g2.chrom == ch).to_numpy()
        if not mk.any(): continue
        x[mk] = g2.mid[mk] + off; centers.append(off + g2.mid[mk].max() / 2); off += g2.mid[mk].max() * 1.02
    g2["x"] = x
    fig, ax = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
    for axi, col, tit in [(ax[0], "p_naive", f"NAIVE (no kinship) λ_GC={lgc_n:.2f}"),
                          (ax[1], "p_kin", f"KINSHIP-CORRECTED (LOCO, common GRM) λ_GC={lgc_k:.2f}  (q<0.05: {nq})")]:
        nlp = -np.log10(np.clip(g2[col], 1e-300, 1))
        for i, ch in enumerate(CHROMS):
            mk = (g2.chrom == ch).to_numpy()
            axi.scatter(g2.x[mk], nlp[mk], s=5, c=["#3b4cc0", "#7aa0c4"][i % 2], alpha=.5, edgecolors="none", rasterized=True)
        axi.axhline(-np.log10(bonf), color="firebrick", lw=.9, ls="--")
        axi.set_ylabel(r"$-\log_{10}p$"); axi.set_title(tit, loc="left", fontsize=10)
        axi.spines[["top", "right"]].set_visible(False)
    ax[1].set_xticks(centers); ax[1].set_xticklabels(CHROMS); ax[1].set_xlabel("genome position")
    fig.suptitle(f"Site {SITE} founder-panel fitness GWAS — all {nF} founders, linear+QN trait, "
                 f"{M:,} haploblocks tested (MAF≥{MAC_MIN/nF:.0%})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig(f"{H}/site{SITE}_founder_gwas231.png", dpi=150)
    print(f"\n[done] {H}/site{SITE}_founder_gwas231.csv + _meta.json + manhattan")


if __name__ == "__main__":
    main()
