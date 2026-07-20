#!/usr/bin/env python
"""Multi-trait founder GWAS across sites — joint / global / climate decomposition (STAGE 1).

Stacks the single-site founder GWAS (founder_gwas_231: linear+QN trait, LOCO GRM from common
markers, EMMAX/P3D) across all sites with persistent gen-1/2/3 plots. Each block then has an
S-vector of per-site effects z_b (genomic-controlled). Multi-trait meta-analysis (Bolormaa 2014)
combines them with the cross-site null covariance C = corr(z over the genome):

  JOINT  (S df) : Q_b = z_b^T C^-1 z_b ~ chi2_S            -- selected at ANY site
  GLOBAL (1 df) : 1^T C^-1 z_b / sqrt(1^T C^-1 1) ~ N(0,1) -- generalist (same direction everywhere)
  CLIMATE(1 df) : c^T C^-1 z_b / sqrt(c^T C^-1 c) ~ N(0,1) -- climate-differential (local adaptation);
                  c = centered site bio1, C-orthogonalized against GLOBAL so the two partition cleanly.

Founder = independent unit (no plot pseudoreplication); site = trait. The rigorous climate null
(site-climate permutation) is STAGE 2; here we report lambda_GC-calibrated results.
Writes multisite_founder_gwas.{csv,npz,_meta.json} + manhattan. Env: kmate.
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

H = "analysis/grenenet_gea/archive/window_hapfreq_retired/hapfreq"
WIN = lib.OUT; SEED = lib.SEEDMIX
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
SUFFIX = os.environ.get("OUT_SUFFIX", "")   # e.g. "_clq90" to keep the K500 and clq0.9 runs side by side
MAC_MIN = int(os.environ.get("MAC_MIN", 3))                # TEST set ~MAF>=1%
MAC_GRM = int(os.environ.get("MAC_GRM", 12))               # GRM set ~MAF>=5% (stable delta)
# Trait = per-founder selection (OLS logit-slope) from p0 over the generations each site HAS, up to
# its last available generation: a single-year site uses p0->gen1, a surviving site uses p0->...->gen3.
# Default "1,2,3" = use the full available trajectory per site (gen1 is required as the anchor).
# CAVEAT: this makes the temporal window vary by site, and window length correlates with climate
# (early-dying sites, mostly hot, only reach gen1 vs gen3 at surviving cold sites) -> a structural
# climate confound the permutation null does NOT remove. Set TRAIT_GENS="1" for the uniform-window
# (gen1-only) sensitivity that is immune to it.
TRAIT_GENS = tuple(int(x) for x in os.environ.get("TRAIT_GENS", "1,2,3").split(","))
# Site 33 is an outlier garden (weird data) -> excluded, matching predict_ecotype_performance.py.
EXCLUDE_SITES = {int(x) for x in os.environ.get("EXCLUDE_SITES", "33").split(",") if x.strip()}


def build_trait_site(sd, p0, gh, nF, use_gens=TRAIT_GENS):
    """Per-founder selection slope from p0 over the requested generations (flower-weighted site-level
    cell per gen, pooled across plots). Needs gen1 (the anchor). Returns (slope, present_gens) or None."""
    cell = {}
    for gen, g in sd.groupby("generation"):
        if int(gen) not in use_gens:
            continue
        hs, ws = [], []
        for _, r in g.iterrows():
            h = gh(str(r.sampleid))
            if h is not None:
                w = r.flowerscollected if np.isfinite(r.flowerscollected) and r.flowerscollected > 0 else 1.0
                hs.append(h); ws.append(w)
        if hs:
            ws = np.asarray(ws); cell[int(gen)] = (np.vstack(hs) * ws[:, None]).sum(0) / ws.sum()
    present = sorted(cell)
    if 1 not in present:                                   # require the gen1 anchor
        return None, []
    t = np.array([0.0] + [float(g) for g in present]); tc = t - t.mean()
    Y = np.vstack([p0] + [np.clip(cell[g], 0.0, 1.0) for g in present])
    slope = (tc[:, None] * Y).sum(0) / (tc @ tc)           # OLS per-generation selection rate
    return slope, present


def loco_emmax(yq, Gp, chrom_m, Gg, chrom_g, nF):
    """Single-trait kinship-corrected EMMAX (LOCO, common GRM). Returns z (M,), delta dict."""
    M = Gp.shape[1]; z = np.full(M, np.nan); dby = {}

    def grm(cols):
        sub = Gg[:, cols]; p = sub.mean(0)
        Z = (sub - p) / np.sqrt(p * (1 - p) + 1e-9)
        return ((Z @ Z.T) / cols.sum() + (Z @ Z.T).T / cols.sum()) / 2

    for ch in CHROMS:
        on = (chrom_m == ch)
        if not on.any():
            continue
        K = grm(chrom_g != ch); lam, U = np.linalg.eigh(K); lam = np.clip(lam, 1e-9, None)
        yr = U.T @ yq; Xr = U.T @ np.ones((nF, 1))
        delta = emma_reml_delta(yr, Xr, lam); dby[ch] = float(delta)
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
        z[on] = b / se
    return z, dby


def bh(pv):
    m = len(pv); o = pv.argsort(); q = np.empty(m)
    q[o] = np.minimum.accumulate((pv[o] * m / (np.arange(m) + 1))[::-1])[::-1]; return np.clip(q, 0, 1)


def lamgc(pv): return np.median(stats.chi2.isf(np.clip(pv, 1e-300, 1), 1)) / stats.chi2.ppf(0.5, 1)


def main():
    G, founders, reg = build_genotype(); nF = len(founders)
    reg["unit"] = reg.chrom + ":" + reg.start.astype(str) + "-" + reg.end.astype(str)
    seeds = sorted({p.split("/")[-1].split("_Chr")[0] for p in glob.glob(f"{SEED}/*_Chr1.h_per_chrom.npz")})
    p0 = np.mean([genome_h(s, SEED) for s in seeds], 0)
    clim = lib.load_climate()
    pt = lib.pool_table()
    pt = pt[pt.sampleid.astype(str).apply(lambda x: os.path.exists(f"{WIN}/{x}_Chr1.h_per_chrom.npz"))]
    gh_cache = {}
    def gh(s):
        if s not in gh_cache: gh_cache[s] = genome_h(s, WIN)
        return gh_cache[s]

    # markers: TEST set MAF>=1%, GRM set common MAF>=5%
    cnt = G.sum(0); blk = reg.block.to_numpy()
    order = np.lexsort((-cnt, blk)); sbk = blk[order]
    is_ref = np.zeros(len(cnt), bool); is_ref[order[np.concatenate([[True], sbk[1:] != sbk[:-1]])]] = True
    poly = (~is_ref) & (cnt >= MAC_MIN) & (cnt <= nF - MAC_MIN)
    grm_set = (~is_ref) & (cnt >= MAC_GRM) & (cnt <= nF - MAC_GRM)
    Gp = G[:, poly].astype(np.float64); regp = reg[poly].reset_index(drop=True); M = Gp.shape[1]
    chrom_m = regp.chrom.to_numpy()
    Gg = G[:, grm_set].astype(np.float64); chrom_g = reg.chrom.to_numpy()[grm_set]
    print(f"{nF} founders | TEST {M:,} markers (MAF>={MAC_MIN/nF:.1%}); GRM {Gg.shape[1]:,} common (MAF>={MAC_GRM/nF:.0%})")

    # ---- per-site single-trait GWAS -> z (genomic-controlled) ----
    print(f"trait generations = {list(TRAIT_GENS)} (gen1 = strongest-selection anchor)")
    Zcols, sites, gens_used, lam_site, traitvec = [], [], [], [], {}
    for site, sd in pt.groupby("site"):
        if int(site) in EXCLUDE_SITES:
            continue
        b1 = float(clim["bio1"].reindex([int(site)]).iloc[0]) if int(site) in clim.index else np.nan
        if not np.isfinite(b1):
            continue
        s_f, present = build_trait_site(sd, p0, gh, nF)
        if s_f is None:
            continue
        yq = stats.norm.ppf((stats.rankdata(s_f) - 0.5) / nF)
        z, _ = loco_emmax(yq, Gp, chrom_m, Gg, chrom_g, nF)
        lg = np.nanmedian(z ** 2) / stats.chi2.ppf(0.5, 1)              # genomic-control factor (diagnostic only)
        # NO genomic control: kinship already corrects structure; stack raw z, record lg to inspect.
        Zcols.append(z); lam_site.append(float(lg))
        sites.append(int(site)); gens_used.append(present); traitvec[int(site)] = s_f
        sg = (s_f - np.median(s_f)) / (stats.median_abs_deviation(s_f) + 1e-12)
        print(f"  site {int(site):>2}: gens {present}, bio1={b1:.1f}C, lambda={lg:.2f}, "
              f"trait max-z={sg.max():.1f} (founder {founders[s_f.argmax()]})")
    Z = np.vstack(Zcols).T                                              # (M, S) genomic-controlled z
    S = Z.shape[1]; sites = np.array(sites); bio1 = clim["bio1"].reindex(sites).to_numpy(float)
    print(f"\nstacked {S} sites x {M:,} blocks; mean per-site lambda {np.mean(lam_site):.2f}")

    # ---- cross-site null covariance C (Bolormaa: from genome-wide z; sparse signal ~ unbiased) ----
    C = np.corrcoef(Z.T); Cinv = np.linalg.pinv(C)
    one = np.ones(S)
    # GLOBAL contrast
    dg = float(one @ Cinv @ one); g = (Z @ Cinv @ one) / np.sqrt(dg)
    # CLIMATE contrast: standardized bio1, C-orthogonalized against GLOBAL (1-vector)
    c0 = (bio1 - bio1.mean()) / bio1.std()
    c = c0 - (float(one @ Cinv @ c0) / dg) * one
    dc = float(c @ Cinv @ c); clim_z = (Z @ Cinv @ c) / np.sqrt(dc)
    # JOINT (omnibus, S df)
    Q = np.einsum("mi,ij,mj->m", Z, Cinv, Z)
    p_joint = stats.chi2.sf(Q, S); p_global = 2 * stats.norm.sf(np.abs(g)); p_clim = 2 * stats.norm.sf(np.abs(clim_z))

    d = regp[["chrom", "start", "end", "unit"]].copy(); d["mac"] = cnt[poly]
    d["chi2_joint"] = Q; d["p_joint"] = p_joint; d["q_joint"] = bh(p_joint)
    d["z_global"] = g; d["p_global"] = p_global; d["q_global"] = bh(p_global)
    d["z_clim"] = clim_z; d["p_clim"] = p_clim; d["q_clim"] = bh(p_clim)
    bonf = 0.05 / M
    d.sort_values("p_joint").to_csv(f"{H}/multisite_founder_gwas{SUFFIX}.csv", index=False)
    np.savez(f"{H}/multisite_founder_gwas{SUFFIX}.npz", Z=Z, sites=sites, bio1=bio1, C=C,
             chrom=regp.chrom.to_numpy(), start=regp.start.to_numpy(), end=regp.end.to_numpy(), mac=cnt[poly])

    rep = {}; com = d.mac >= MAC_GRM                       # common-marker mask (rare-artifact guard)
    for tag, pc, qc in [("JOINT ", "p_joint", "q_joint"), ("GLOBAL", "p_global", "q_global"), ("CLIMATE", "p_clim", "q_clim")]:
        nb = int((d[pc] < bonf).sum()); nq = int((d[qc] < 0.05).sum())
        nqc = int((d[qc] < 0.05).loc[com].sum())           # FDR hits among common (MAF>=5%) markers
        rep[tag.strip()] = dict(lam=float(lamgc(d[pc].values)), n_bonf=nb, n_fdr=nq, n_fdr_common=nqc, best_q=float(d[qc].min()))
        print(f"  [{tag}] lambda_GC={lamgc(d[pc].values):.2f} | Bonferroni {nb} | FDR q<0.05 {nq} "
              f"({nqc} common MAC>={MAC_GRM}) | best q={d[qc].min():.3g}")

    meta = dict(n_founders=nF, M=int(M), S=int(S), sites=sites.tolist(), excluded_sites=sorted(EXCLUDE_SITES),
                memb_tag=os.environ.get("MEMB_TAG", "K500"),
                trait_gens=list(TRAIT_GENS),
                gens_used=gens_used, bio1=bio1.tolist(), per_site_lambda=lam_site, mac_min=MAC_MIN, mac_grm=MAC_GRM,
                mean_cross_site_C=float(C[np.triu_indices(S, 1)].mean()), tests=rep,
                model="multi-trait founder GWAS (Bolormaa meta over single-site EMMAX): joint/global/climate; STAGE1 lambda-calibrated, perm null TODO")
    json.dump(meta, open(f"{H}/multisite_founder_gwas{SUFFIX}_meta.json", "w"), indent=2)

    # ---- Manhattan: joint / global / climate ----
    g2 = d.copy(); g2["mid"] = (g2.start + g2.end) / 2
    g2 = g2.sort_values(["chrom", "start"]).reset_index(drop=True)
    off, centers, x = 0.0, [], np.zeros(len(g2))
    for ch in CHROMS:
        mk = (g2.chrom == ch).to_numpy()
        if not mk.any(): continue
        x[mk] = g2.mid[mk] + off; centers.append(off + g2.mid[mk].max() / 2); off += g2.mid[mk].max() * 1.02
    g2["x"] = x
    fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    for axi, col, tit in [(ax[0], "p_joint", f"JOINT (any-site selection, {S}df)  λ={rep['JOINT']['lam']:.2f}  q<.05: {rep['JOINT']['n_fdr']}"),
                          (ax[1], "p_global", f"GLOBAL (generalist)  λ={rep['GLOBAL']['lam']:.2f}  q<.05: {rep['GLOBAL']['n_fdr']}"),
                          (ax[2], "p_clim", f"CLIMATE (local adaptation, bio1)  λ={rep['CLIMATE']['lam']:.2f}  q<.05: {rep['CLIMATE']['n_fdr']}")]:
        nlp = -np.log10(np.clip(g2[col], 1e-300, 1))
        for i, ch in enumerate(CHROMS):
            mk = (g2.chrom == ch).to_numpy()
            axi.scatter(g2.x[mk], nlp[mk], s=5, c=["#3b4cc0", "#7aa0c4"][i % 2], alpha=.5, edgecolors="none", rasterized=True)
        axi.axhline(-np.log10(bonf), color="firebrick", lw=.9, ls="--")
        axi.set_ylabel(r"$-\log_{10}p$"); axi.set_title(tit, loc="left", fontsize=10); axi.spines[["top", "right"]].set_visible(False)
    ax[2].set_xticks(centers); ax[2].set_xticklabels(CHROMS); ax[2].set_xlabel("genome position")
    fig.suptitle(f"Multi-trait founder GWAS across {S} sites — joint / global / climate decomposition", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig(f"{H}/multisite_founder_gwas{SUFFIX}.png", dpi=150)
    print(f"\n[done] {H}/multisite_founder_gwas{SUFFIX}.csv + .npz + _meta.json + manhattan")


if __name__ == "__main__":
    main()
