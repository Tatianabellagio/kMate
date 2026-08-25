#!/usr/bin/env python
"""Per-site class-split founder GWAS run with GEMMA instead of the in-house EMMAX/P3D scan.

Decision (user, 2026-08-25): use the community-validated tool rather than our own association
model, so the method is citable in one sentence. This replaces the PER-SITE scan of
class_split_gwas.py with GEMMA 0.98.5 (`-lmm 1`, Wald). Everything AROUND the per-site scan is
unchanged and deliberately so:
  * the trait (selection_s_matrix.npz, rank-normalized per site),
  * the pooled per-chrom LOCO GRM (MAC>=12, call>=90%, ALL classes pooled -- shared across
    classes so a class difference is attributable to the test markers, not the correction),
  * the test-marker filters (MAC>=5, call>=90%) and the snp/nonsnp/sv split,
  * the cross-site meta (Bolormaa et al. 2014 JOINT/GLOBAL/CLIMATE), which GEMMA cannot do:
    GCTA/EMMAX have no multi-phenotype association test and GEMMA's mvLMM cannot identify 30
    traits from 231 founders (930 variance parameters). The meta is a published FORMULA, so it
    needs a citation, not a tool.
GEMMA validation of the estimator being replaced: pearson(z)=0.9975 (gemma_validation.py).

Output schema is a DROP-IN match for class_gwas_{snp,nonsnp,sv}.npz so every downstream consumer
(block collapse, gene attribution, GO enrichment, both notebooks, plot_class_gwas_pngs.py) works
unchanged -- written with a _gemma suffix so the in-house arm stays on disk for comparison.

Two-stage, because a single process would serialize ~3,400 GEMMA runs:
  stage `scan`  (sbatch array, one task per chrom x class = 15 tasks): writes the BIMBAM
                genotype ONCE per chrom x class and reuses it across all 30 sites (the phenotype
                is only a `-n` column switch), then runs GEMMA per site x chunk ->
                results/gemma_gwas/parts/{chrom}_{class}.npz  (pos, Z[M,S], lam_persite)
  stage `meta`  (single task, after the array): concatenates chroms, applies the same
                finite-marker filter and Bolormaa meta as class_split_gwas.py, and writes
                results/gemma_gwas/class_gwas_{class}_gemma.npz + class_gwas_summary_gemma.json

Usage:
    python gemma_persite_gwas.py scan --chrom chr1 --class snp
    python gemma_persite_gwas.py meta
Env: `kmate` python + the `gemma` binary from env gwas_tools on PATH. See run_gemma_gwas.sbatch.
"""
from __future__ import annotations
import os, sys, json, time, shutil, argparse, subprocess
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from class_split_gwas import (_panel_order, _load_chrom_212, _classes, _dense_imputed, _ZZt,
                              bolormaa, CHROMS, MIN_MAC_TEST, MIN_MAC_GRM, CALL_MIN, WIN)
from gemma_validation import GEMMA_MAX, build_pooled_grms

OUT = f"{lib.GEA}/r3_persite_gwas/results/gemma_gwas"
PARTS = f"{OUT}/parts"
IO = f"{OUT}/gemma_io"
CLASSES = ("snp", "nonsnp", "sv")


def _pheno_path(Smat, N):
    """Rank-normalized per-site phenotype matrix, one column per site (GEMMA -n is 1-based)."""
    os.makedirs(IO, exist_ok=True)
    path = f"{IO}/pheno_all.txt"
    if not os.path.exists(path):
        P = np.column_stack([stats.norm.ppf((stats.rankdata(Smat[si]) - 0.5) / N)
                             for si in range(Smat.shape[0])])
        np.savetxt(path, P, fmt="%.10f")
    return path


def scan(chrom, clsname):
    t0 = time.time()
    if shutil.which("gemma") is None:
        sys.exit("gemma not on PATH -- prepend $HOME/miniforge3/envs/gwas_tools/bin")
    os.makedirs(PARTS, exist_ok=True); os.makedirs(IO, exist_ok=True)
    order, Smat, sites, bio1 = _panel_order()
    N, S = len(order), len(sites)

    Kchrom, Ksum, Mtot = build_pooled_grms(order, N)
    K_loco = (Ksum - Kchrom[chrom][0]) / max(Mtot - Kchrom[chrom][1], 1)
    kin = f"{IO}/kin_{chrom}.txt"
    if not os.path.exists(kin):
        np.savetxt(kin, K_loco, fmt="%.10f")
    pheno = _pheno_path(Smat, N)

    pos, rl, al, vp, vc, n_alt, n_cal = _load_chrom_212(chrom, order)
    mac = np.minimum(n_alt, N - n_alt)
    vclass = _classes(rl, al)
    keep = (mac >= MIN_MAC_TEST) & (n_cal >= CALL_MIN * N)
    want = {"snp": vclass == 0, "nonsnp": vclass != 0, "sv": vclass == 2}[clsname]
    idx = np.where(keep & want)[0]
    if len(idx) == 0:
        np.savez(f"{PARTS}/{chrom}_{clsname}.npz", pos=np.array([]), Z=np.zeros((0, S)),
                 lam_persite=np.zeros(S)); return
    marker_pos = pos[idx]
    print(f"{chrom}/{clsname}: {len(idx):,} test markers, {S} sites, {N} founders", flush=True)

    # BIMBAM written ONCE, reused for every site. rs carries the ordinal -- NOT chrom:pos, which
    # would be ambiguous at the 2.14% of positions holding >1 biallelic record (multiallelic
    # decomposition); see gemma_validation.py.
    chunks = np.array_split(np.arange(len(idx)), int(np.ceil(len(idx) / GEMMA_MAX)))
    for ci, ch in enumerate(chunks):
        gp = f"{IO}/geno_{chrom}_{clsname}_{ci}.bimbam"
        if os.path.exists(gp):
            continue
        g, _ = _dense_imputed(vp, vc, n_alt, n_cal, idx[ch])
        gt = g.T
        with open(gp, "w") as fh:
            for k, j in enumerate(ch):
                fh.write(f"{chrom}:{marker_pos[j]}:{j},A,T," +
                         ",".join(f"{v:.6g}" for v in gt[k]) + "\n")
        del g, gt
    print(f"  {len(chunks)} chunk(s) <= {GEMMA_MAX:,} markers ({time.time()-t0:.0f}s)", flush=True)

    Z = np.full((len(idx), S), np.nan)
    lam = np.zeros(S)
    for si in range(S):
        zs = np.full(len(idx), np.nan)
        for ci in range(len(chunks)):
            pref = f"{chrom}_{clsname}_s{si}_c{ci}"
            cmd = ["gemma", "-g", f"{IO}/geno_{chrom}_{clsname}_{ci}.bimbam", "-p", pheno,
                   "-k", kin, "-lmm", "1", "-n", str(si + 1),
                   # QC wide open: GEMMA must test EXACTLY our marker set, or a marker-count
                   # difference confounds the comparison with the in-house arm. GEMMA DOES
                   # filter by default (-maf 0.01 -miss 0.05 -r2 0.9999) -- we disable it
                   # because ALL QC is already done upstream, and GEMMA's own MAF is WRONG for
                   # this data (measured 2026-08-25):
                   #   * GEMMA assumes diploid dosage, so it computes maf = mean(geno)/2 = p/2
                   #     for our haploid 0/1 coding. Verified: predicted vs actual drops at
                   #     -maf 0.02 = 2921 vs 2927, at -maf 0.05 = 5635 vs 5634.
                   #   * That makes it ONE-SIDED -- it can never exceed 0.5, so it cannot catch
                   #     a near-fixed marker. A probe with 230/231 founders ALT (true MAC=1)
                   #     PASSED GEMMA's default QC, read as maf 0.4978.
                   #   * A fully monomorphic probe was caught, but by -r2 0.9999 (constant
                   #     vector collinear with the intercept), NOT by -maf -- and -r2 1 turns
                   #     that off too.
                   # Our upstream filter is symmetric and strictly stronger: MAC =
                   # min(n_alt, N-n_alt) >= 5 (MAF 2.16%) excludes monomorphic (MAC 0) and
                   # near-fixed (MAC 1) before GEMMA sees them. Nothing degenerate can reach it.
                   # (At GEMMA's DEFAULT -maf 0.01 it would drop 0 of our markers anyway -- our
                   # MAC>=5 floor lands at maf 0.0108 -- but that is a 0.0008 margin, i.e. a
                   # coincidence, so we do not depend on it.)
                   "-maf", "0", "-miss", "1", "-r2", "1", "-hwe", "0",
                   "-outdir", IO, "-o", pref]
            r = subprocess.run(cmd, capture_output=True, text=True)
            assoc = f"{IO}/{pref}.assoc.txt"
            if not os.path.exists(assoc):
                print(r.stdout[-1500:]); print(r.stderr[-1500:])
                sys.exit(f"GEMMA failed: {pref}")
            a = np.genfromtxt(assoc, names=True, dtype=None, encoding=None)
            gj = np.array([int(str(s).split(":")[2]) for s in np.atleast_1d(a["rs"])])
            zs[gj] = np.atleast_1d(a["beta"]) / np.atleast_1d(a["se"])
            os.remove(assoc); os.remove(f"{IO}/{pref}.log.txt")
        Z[:, si] = zs
        lam[si] = np.nanmedian(zs ** 2) / stats.chi2.ppf(0.5, 1)
        print(f"  site {sites[si]} ({si+1}/{S}): lambda={lam[si]:.3f} ({time.time()-t0:.0f}s)",
              flush=True)
    np.savez(f"{PARTS}/{chrom}_{clsname}.npz", pos=marker_pos, Z=Z, lam_persite=lam)
    for ci in range(len(chunks)):
        os.remove(f"{IO}/geno_{chrom}_{clsname}_{ci}.bimbam")     # large; parts npz is the keep
    print(f"wrote {PARTS}/{chrom}_{clsname}.npz ({time.time()-t0:.0f}s)", flush=True)


def meta():
    """Identical cross-site meta to class_split_gwas.py, on the GEMMA per-site Z."""
    order, Smat, sites, bio1 = _panel_order()
    summary = {"n_founders": len(order), "n_sites": len(sites), "sites": sites.tolist(),
               "estimator": "gemma-0.98.5 -lmm 1 (Wald), per-site; meta = Bolormaa et al. 2014"}
    win_stat = {}
    for clsname in CLASSES:
        pos_all, chrom_all, Zs, lams = [], [], [], []
        for cl in CHROMS:
            p = f"{PARTS}/{cl}_{clsname}.npz"
            if not os.path.exists(p):
                sys.exit(f"missing {p} -- run the scan stage for {cl}/{clsname} first")
            z = np.load(p, allow_pickle=True)
            if len(z["pos"]) == 0:
                continue
            pos_all.append(z["pos"]); chrom_all.append(np.array([cl] * len(z["pos"])))
            Zs.append(z["Z"]); lams.append(z["lam_persite"])
        pos_all = np.concatenate(pos_all); chrom_all = np.concatenate(chrom_all)
        Zall = np.concatenate(Zs, axis=0)
        finite = np.isfinite(Zall).all(1)
        if (~finite).sum():
            print(f"[{clsname}] dropping {int((~finite).sum())} markers not testable at every "
                  f"site before the meta", flush=True)
        pos_all, chrom_all, Zall = pos_all[finite], chrom_all[finite], Zall[finite]
        Q, p_joint, z_global, p_global, z_clim, p_clim, C = bolormaa(Zall, sites, bio1)
        q_joint, q_global, q_clim = lib.bh(p_joint), lib.bh(p_global), lib.bh(p_clim)
        M = len(pos_all); bonf = 0.05 / M
        np.savez(f"{OUT}/class_gwas_{clsname}_gemma.npz", chrom=chrom_all, pos=pos_all, Z=Zall,
                 chi2_joint=Q, p_joint=p_joint, q_joint=q_joint,
                 z_global=z_global, p_global=p_global, q_global=q_global,
                 z_clim=z_clim, p_clim=p_clim, q_clim=q_clim,
                 lam_persite=np.mean(lams, axis=0), sites=sites, bio1=bio1)
        win = (pos_all // WIN).astype(np.int64)
        key = np.array([f"{c}:{w}" for c, w in zip(chrom_all, win)])
        nlp = -np.log10(np.clip(p_joint, 1e-300, 1))
        winmax = {}
        for k, v in zip(key, nlp):
            if k not in winmax or v > winmax[k]:
                winmax[k] = v
        win_stat[clsname] = winmax
        summary[clsname] = dict(
            n_markers=M, lambda_JOINT=lib.lamgc(p_joint), lambda_GLOBAL=lib.lamgc(p_global),
            lambda_CLIM=lib.lamgc(p_clim), mean_persite_lambda=float(np.mean(lams)),
            bonferroni_joint=int((p_joint < bonf).sum()), fdr_joint=int((q_joint < 0.05).sum()),
            fdr_global=int((q_global < 0.05).sum()), fdr_clim=int((q_clim < 0.05).sum()),
            best_p_joint=float(p_joint.min()))
        print(f"[{clsname}] M={M:,} lambda_JOINT={summary[clsname]['lambda_JOINT']:.3f} "
              f"Bonf={summary[clsname]['bonferroni_joint']} "
              f"FDR={summary[clsname]['fdr_joint']} best_p={p_joint.min():.2e}", flush=True)

    def compare(c1, c2):
        shared = sorted(set(win_stat[c1]) & set(win_stat[c2]))
        a = np.array([win_stat[c1][k] for k in shared]); b = np.array([win_stat[c2][k] for k in shared])
        rho, rp = stats.spearmanr(a, b)
        out = {"shared_windows_n": len(shared), "window_spearman_rho": float(rho),
               "window_spearman_p": float(rp)}
        for topn in (0.001, 0.005, 0.01):
            k1 = set(sorted(win_stat[c1], key=lambda k: -win_stat[c1][k])[:max(1, int(topn * len(win_stat[c1])))])
            k2 = set(sorted(win_stat[c2], key=lambda k: -win_stat[c2][k])[:max(1, int(topn * len(win_stat[c2])))])
            out[f"top{topn}_jaccard"] = len(k1 & k2) / max(len(k1 | k2), 1)
        return out

    summary["comparisons"] = {f"snp_vs_{c}": compare("snp", c) for c in ("nonsnp", "sv")}
    summary.update(summary["comparisons"]["snp_vs_nonsnp"])
    json.dump(summary, open(f"{OUT}/class_gwas_summary_gemma.json", "w"), indent=2, default=str)
    print(f"wrote {OUT}/class_gwas_{{snp,nonsnp,sv}}_gemma.npz + class_gwas_summary_gemma.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["scan", "meta"])
    ap.add_argument("--chrom", choices=CHROMS)
    ap.add_argument("--class", dest="cls", choices=CLASSES)
    a = ap.parse_args()
    if a.stage == "scan":
        if not a.chrom or not a.cls:
            sys.exit("scan needs --chrom and --class")
        scan(a.chrom, a.cls)
    else:
        meta()
