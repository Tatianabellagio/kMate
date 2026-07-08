#!/usr/bin/env python
"""Split single-marker founder GWAS by class (SNP vs non-SNP) -- do SNP-based and non-SNP-based
association scans find the same PEAKS, and how does variance explained compare, per-site and in
the cross-site JOINT/GLOBAL/CLIMATE meta?

The variance-partition work (varexp_selection.py / varexp_untagged.py) answers a KINSHIP-level
question: do SNP and non-SNP GRMs explain the same total selection-trait variance (yes, even for
SNP-untagged non-SNP markers -- both just reconstruct the same founder genealogy). This script
answers the complementary PER-MARKER question at raw-variant resolution (not the coarser
clq0.9-block/hap-cluster unit `founder_gwas_multisite.py` uses): if you run an actual
single-marker association scan with only SNPs vs only non-SNP (indel+SV) markers, do you get the
same top hits, and the same genomic-control-calibrated inflation (a per-scan "how much signal did
this marker class pick up" proxy)?

Method (mirrors founder_gwas_multisite.py's loco_emmax + Bolormaa meta, adapted to be class-split
and to run efficiently at ~2.3M-marker scale):
  - Trait: the LOCKED selection coefficient `s` (selection_s_matrix.npz, 31 sites x 212
    analyzable founders) -- same trait as the variance-partition work, rank-normalized per site.
  - GRM: ONE shared LOCO kinship correction per chromosome, built from ALL classes pooled
    (MAC>=12, call>=90%, matches build_class_grms.py / multisite script's "common marker" floor)
    -- since K_snp/K_nonsnp are near-identical (corr up to 0.998, see varexp thread), a shared
    correction isolates the TEST-marker-class effect instead of also varying the correction.
    LOCO done via the standard total-minus-chromosome trick (no repeated whole-genome loads).
  - Test markers: MAC>=5, call>=90% (matches class_grms.py), split SNP vs non-SNP (indel+SV
    pooled, the primary 2-way split already locked with the user) vs strict SV-only
    (|alt_len-ref_len|>50bp, a subset of non-SNP -- the secondary 3-way breakdown).
  - Per (chrom, class): ONE eigendecomposition + ONE U.T-rotation of the test-marker genotype
    block, reused across all 31 sites (only the REML delta + GLS step re-runs per site) --
    this is what keeps ~2.3M markers x 31 sites tractable interactively.
  - Per-site z genomic-control calibrated (z / sqrt(lambda_GC)), then stacked into per-class
    Z[M,S] and run through the identical Bolormaa JOINT (S df, any-site) / GLOBAL (1df,
    generalist) / CLIMATE (1df, bio1-differential) formulas as founder_gwas_multisite.py.
  - Peaks comparison: Bonferroni/FDR hit counts + lambda_GC per class; genome binned into 20kb
    windows, per-window max -log10(p_joint) correlated between classes (Spearman) as the
    "same regions light up" test; top-N hit window overlap (Jaccard).

Reads results/grenenet_gea/varexp/selection_s_matrix.npz + panel/arch3 var_pa per chrom.
Writes results/grenenet_gea/varexp/class_gwas_{snp,nonsnp}.npz + class_gwas_summary.json.
Env: kmate. Compute-node only (dense per-chrom marker blocks).
"""
from __future__ import annotations
import os, sys, json, time, argparse
import numpy as np
import scipy.sparse as sp
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from founder_genotype import emma_reml_delta

OUT = f"{lib.GEA}/varexp"
CHROMS = ["chr1", "chr2", "chr3", "chr4", "chr5"]
MIN_MAC_TEST = 5
MIN_MAC_GRM = 12
CALL_MIN = 0.9
WIN = 20_000          # window size (bp) for the peak-overlap comparison


def _panel_order():
    S = np.load(f"{OUT}/selection_s_matrix.npz", allow_pickle=True)
    tf = S["founders"].astype(str)
    ana = S["analyzable"].astype(bool)
    panel = np.load(f"{lib.PROJ}/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz",
                    allow_pickle=True)["founders"].astype(str)
    gpos = {f: i for i, f in enumerate(panel)}
    keep = np.array([ana[i] and tf[i] in gpos for i in range(len(tf))])
    order = np.array([gpos[tf[i]] for i in range(len(tf)) if keep[i]])
    return order, S["S"][:, keep], S["sites"], S["bio1"]


def _load_chrom_212(cl, order):
    base = f"{lib.PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
    meta = np.load(f"{base}.meta.npz", allow_pickle=True)
    vp = sp.load_npz(f"{base}.var_pa.npz").tocsr()[order].tocsc()
    vc = sp.load_npz(f"{base}.var_called.npz").tocsr()[order].tocsc()
    n_alt = np.asarray(vp.sum(0)).ravel().astype(float)
    n_cal = np.asarray(vc.sum(0)).ravel().astype(float)
    return (meta["pos"].astype(np.int64), meta["ref_len"].astype(int), meta["alt_len"].astype(int),
            vp, vc, n_alt, n_cal)


def _classes(rl, al):
    dlen = np.abs(al - rl)
    return np.where((rl == 1) & (al == 1), 0, np.where(dlen > 50, 2, 1))


def _dense_imputed(vp, vc, n_alt, n_cal, cols):
    g = vp[:, cols].toarray().astype(np.float64)
    called = vc[:, cols].toarray().astype(bool)
    p = np.where(n_cal[cols] > 0, n_alt[cols] / np.maximum(n_cal[cols], 1), 0.0)
    miss = ~called
    if miss.any():
        g[miss] = np.take(p, np.where(miss)[1])
    return g, p


def _ZZt(g, p):
    sd = np.sqrt(np.maximum(p * (1 - p), 1e-6))
    Z = (g - p) / sd
    return Z @ Z.T


def bolormaa(Z, sites, bio1):
    S = Z.shape[1]
    C = np.corrcoef(Z.T); Cinv = np.linalg.pinv(C)
    one = np.ones(S)
    dg = float(one @ Cinv @ one)
    z_global = (Z @ Cinv @ one) / np.sqrt(dg)
    c0 = (bio1 - bio1.mean()) / bio1.std()
    c = c0 - (float(one @ Cinv @ c0) / dg) * one
    dc = float(c @ Cinv @ c)
    z_clim = (Z @ Cinv @ c) / np.sqrt(dc)
    Q = np.einsum("mi,ij,mj->m", Z, Cinv, Z)
    p_joint = stats.chi2.sf(Q, S)
    p_global = 2 * stats.norm.sf(np.abs(z_global))
    p_clim = 2 * stats.norm.sf(np.abs(z_clim))
    return Q, p_joint, z_global, p_global, z_clim, p_clim, C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true", help="chr1 only, 3 sites -- fast correctness check")
    args = ap.parse_args()
    t0 = time.time()

    order, Smat, sites, bio1 = _panel_order()
    N = len(order)
    chroms = ["chr1"] if args.smoke else CHROMS
    site_idx = range(3) if args.smoke else range(len(sites))
    sites_use = sites[list(site_idx)]
    bio1_use = bio1[list(site_idx)]
    print(f"{N} analyzable founders, {len(sites_use)} sites{' [SMOKE]' if args.smoke else ''}", flush=True)

    # ---- pass 1: per-chrom GRM partial sums (pooled all classes) + test-marker class masks ----
    chrom_data = {}
    Ksum_total = np.zeros((N, N)); Mtot = 0
    Kchrom = {}
    for cl in chroms:
        pos, rl, al, vp, vc, n_alt, n_cal = _load_chrom_212(cl, order)
        mac = np.minimum(n_alt, N - n_alt)
        vclass = _classes(rl, al)
        keep_grm = (mac >= MIN_MAC_GRM) & (n_cal >= CALL_MIN * N)
        gcols = np.where(keep_grm)[0]
        g, p = _dense_imputed(vp, vc, n_alt, n_cal, gcols)
        K = _ZZt(g, p); m = len(gcols)
        Kchrom[cl] = (K, m)
        Ksum_total += K; Mtot += m

        keep_test = (mac >= MIN_MAC_TEST) & (n_cal >= CALL_MIN * N)
        snp_idx = np.where(keep_test & (vclass == 0))[0]
        nonsnp_idx = np.where(keep_test & (vclass != 0))[0]
        sv_idx = np.where(keep_test & (vclass == 2))[0]      # strict SV only (|dlen|>50bp)
        chrom_data[cl] = dict(pos=pos, vp=vp, vc=vc, n_alt=n_alt, n_cal=n_cal,
                              snp_idx=snp_idx, nonsnp_idx=nonsnp_idx, sv_idx=sv_idx)
        print(f"{cl}: grm={m:,} snp_test={len(snp_idx):,} nonsnp_test={len(nonsnp_idx):,} "
              f"sv_test={len(sv_idx):,} ({time.time()-t0:.0f}s)", flush=True)
    print(f"genome GRM markers total = {Mtot:,}", flush=True)

    # ---- pass 2: LOCO GWAS per (chrom, class); one eigh + one rotation reused across all sites ----
    CLASSES = ("snp", "nonsnp", "sv")
    per_class = {c: [] for c in CLASSES}
    lam_track = {c: [] for c in CLASSES}
    for cl in chroms:
        d = chrom_data[cl]
        denom = max(Mtot - Kchrom[cl][1], 1)
        K_loco = (Ksum_total - Kchrom[cl][0]) / denom
        lam_eig, U = np.linalg.eigh(K_loco); lam_eig = np.clip(lam_eig, 1e-9, None)
        Xr = (U.T @ np.ones((N, 1)))[:, 0]
        for clsname, idx in (("snp", d["snp_idx"]), ("nonsnp", d["nonsnp_idx"]), ("sv", d["sv_idx"])):
            if len(idx) == 0:
                continue
            g, p = _dense_imputed(d["vp"], d["vc"], d["n_alt"], d["n_cal"], idx)
            Br = U.T @ g                                    # [N x n_markers] -- reused across sites
            Zmat = np.full((len(idx), len(sites_use)), np.nan)
            lam_site = np.zeros(len(sites_use))
            for si in range(len(sites_use)):
                y = Smat[list(site_idx)[si]]
                yq = stats.norm.ppf((stats.rankdata(y) - 0.5) / N)
                yr = U.T @ yq
                delta = emma_reml_delta(yr, Xr[:, None], lam_eig)
                w = 1.0 / (lam_eig + delta)
                Saa = (w * Xr * Xr).sum(); Say = (w * Xr * yr).sum()
                Sab = ((w * Xr)[None, :] @ Br).ravel()
                Sbb = (w[:, None] * Br ** 2).sum(0)
                Sby = ((w * yr)[None, :] @ Br).ravel()
                det = Saa * Sbb - Sab ** 2
                # near-singular guard: a handful of markers per (chrom,site) have a rotated
                # genotype vector Br numerically collinear with the rotated intercept Xr
                # (relative det ~1e-16, machine-epsilon degenerate -- vs >1e-6 for every
                # legitimate marker, a clean 10-order-of-magnitude gap). Uncaught, these blow
                # up to |z|~1e19 and single-handedly wreck the cross-site corr(Z) used by the
                # Bolormaa meta for ALL markers. Mask them out as untestable (NaN).
                singular = det <= 1e-8 * np.maximum(Saa * Sbb, 1e-300)
                detc = np.where(singular, np.nan, det)
                b = (Saa * Sby - Sab * Say) / detc
                alpha = (Sbb * Say - Sab * Sby) / detc
                rss = (w * yr * yr).sum() - alpha * Say - b * Sby
                s2 = np.clip(rss, 1e-30, None) / (N - 2)
                se = np.sqrt(np.clip(s2 * Saa / detc, 1e-30, None))
                z = b / se
                lg = np.nanmedian(z ** 2) / stats.chi2.ppf(0.5, 1)
                lam_site[si] = lg
                # NO genomic control: the LOCO kinship GRM already corrects structure. Record lg
                # as a diagnostic only and stack the RAW kinship-corrected z. (Per-site GC was
                # inconsistent with founder_gwas_231, and self-inflated wherever lambda<1.)
                Zmat[:, si] = z
            per_class[clsname].append((d["pos"][idx], np.array([cl] * len(idx)), Zmat))
            lam_track[clsname].append(lam_site)
            print(f"  {cl} {clsname}: {len(idx):,} markers done ({time.time()-t0:.0f}s)", flush=True)

    # ---- assemble + Bolormaa meta per class ----
    summary = {"n_founders": N, "n_sites": len(sites_use), "sites": sites_use.tolist()}
    win_stat = {}
    for clsname in CLASSES:
        pos_all = np.concatenate([x[0] for x in per_class[clsname]])
        chrom_all = np.concatenate([x[1] for x in per_class[clsname]])
        Zall = np.concatenate([x[2] for x in per_class[clsname]], axis=0)
        finite = np.isfinite(Zall).all(1)
        n_dropped = int((~finite).sum())
        if n_dropped:
            print(f"[{clsname}] dropping {n_dropped} markers with a near-singular GLS fit "
                  f"(untestable) before the cross-site meta", flush=True)
        pos_all, chrom_all, Zall = pos_all[finite], chrom_all[finite], Zall[finite]
        Q, p_joint, z_global, p_global, z_clim, p_clim, C = bolormaa(Zall, sites_use, bio1_use)
        q_joint = lib.bh(p_joint); q_global = lib.bh(p_global); q_clim = lib.bh(p_clim)
        M = len(pos_all); bonf = 0.05 / M
        np.savez(f"{OUT}/class_gwas_{clsname}.npz", chrom=chrom_all, pos=pos_all, Z=Zall,
                 chi2_joint=Q, p_joint=p_joint, q_joint=q_joint,
                 z_global=z_global, p_global=p_global, q_global=q_global,
                 z_clim=z_clim, p_clim=p_clim, q_clim=q_clim,
                 lam_persite=np.mean(lam_track[clsname], axis=0), sites=sites_use, bio1=bio1_use)
        win = (pos_all // WIN).astype(np.int64)
        key = np.array([f"{c}:{w}" for c, w in zip(chrom_all, win)])
        nlp = -np.log10(np.clip(p_joint, 1e-300, 1))
        order_nlp = np.argsort(-nlp)
        winmax = {}
        for k, v in zip(key, nlp):
            if k not in winmax or v > winmax[k]:
                winmax[k] = v
        win_stat[clsname] = winmax
        summary[clsname] = dict(
            n_markers=M, lambda_JOINT=lib.lamgc(p_joint), lambda_GLOBAL=lib.lamgc(p_global),
            lambda_CLIM=lib.lamgc(p_clim), mean_persite_lambda=float(np.mean(lam_track[clsname])),
            bonferroni_joint=int((p_joint < bonf).sum()), fdr_joint=int((q_joint < 0.05).sum()),
            fdr_global=int((q_global < 0.05).sum()), fdr_clim=int((q_clim < 0.05).sum()),
            best_p_joint=float(p_joint.min()),
            top10_windows=[key[i] for i in order_nlp[:10]])
        print(f"[{clsname}] M={M:,} lambda_JOINT={summary[clsname]['lambda_JOINT']:.3f} "
              f"Bonf={summary[clsname]['bonferroni_joint']} FDR={summary[clsname]['fdr_joint']} "
              f"best_p={summary[clsname]['best_p_joint']:.2e}", flush=True)

    # ---- peaks comparison: shared-window correlation + top-hit overlap, for each pair vs SNP ----
    def compare(c1, c2):
        shared = sorted(set(win_stat[c1]) & set(win_stat[c2]))
        a = np.array([win_stat[c1][k] for k in shared])
        b = np.array([win_stat[c2][k] for k in shared])
        rho, rho_p = stats.spearmanr(a, b)
        out = {"shared_windows_n": len(shared), "window_spearman_rho": float(rho),
               "window_spearman_p": float(rho_p)}
        for topn in (0.001, 0.005, 0.01):
            k1 = set(sorted(win_stat[c1], key=lambda k: -win_stat[c1][k])[:max(1, int(topn * len(win_stat[c1])))])
            k2 = set(sorted(win_stat[c2], key=lambda k: -win_stat[c2][k])[:max(1, int(topn * len(win_stat[c2])))])
            out[f"top{topn}_jaccard"] = len(k1 & k2) / max(len(k1 | k2), 1)
        return out

    comparisons = {}
    for c2 in ("nonsnp", "sv"):
        cmp = compare("snp", c2)
        comparisons[f"snp_vs_{c2}"] = cmp
        print(f"\n[snp vs {c2}] shared windows={cmp['shared_windows_n']}  "
              f"spearman={cmp['window_spearman_rho']:.3f} (p={cmp['window_spearman_p']:.2g})")
        for topn in (0.001, 0.005, 0.01):
            print(f"  top-{topn:.1%} window Jaccard = {cmp[f'top{topn}_jaccard']:.3f}")
    summary["comparisons"] = comparisons
    # back-compat top-level keys (existing notebook reads these for the snp-vs-nonsnp pair)
    summary.update({k: v for k, v in comparisons["snp_vs_nonsnp"].items()})
    json.dump(summary, open(f"{OUT}/class_gwas_summary.json", "w"), indent=2, default=str)
    print(f"\nwrote {OUT}/class_gwas_{{snp,nonsnp,sv}}.npz + class_gwas_summary.json ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
