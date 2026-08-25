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
  * the test-marker filters (MAC>=5, call>=90%) and the snp/nonsnp/sv split.
GEMMA validation of the estimator being replaced: pearson(z)=0.9975 (gemma_validation.py).

**SCOPE (user decision, 2026-08-25): PER-SITE ONLY -- there is no cross-site meta.** The
all-sites question is answered by the GEA track, so this analysis is 30 independent per-garden
LMMs and nothing else. The previous Bolormaa et al. 2014 JOINT/GLOBAL/CLIMATE meta is dropped
rather than ported. That deletes the only stage of this pipeline that was not a named,
community-validated tool: GCTA has no multi-phenotype association test (`--reml-bivar` is
exactly 2 traits; `--mtcojo` conditions summary stats and needs an external LD reference) and
GEMMA's mvLMM cannot identify 30 traits from 231 founders. Every number this pipeline now
produces comes out of GEMMA. The in-house meta lives on in class_split_gwas.py for the record.

Output schema is a DROP-IN match for class_gwas_{snp,nonsnp,sv}.npz so every downstream consumer
(block collapse, gene attribution, GO enrichment, both notebooks, plot_class_gwas_pngs.py) works
unchanged -- written with a _gemma suffix so the in-house arm stays on disk for comparison.

Two-stage, because a single process would serialize ~3,400 GEMMA runs:
  stage `scan`  (sbatch array, one task per chrom x class = 15 tasks): writes the BIMBAM
                genotype ONCE per chrom x class and reuses it across all 30 sites (the phenotype
                is only a `-n` column switch), then runs GEMMA per site x chunk ->
                results/gemma_gwas/parts/{chrom}_{class}.npz  (pos, Z[M,S], lam_persite)
  stage `assemble` (single task, after the array): concatenates chroms per class and attaches
                per-GARDEN significance (two-sided p from the GEMMA Wald z, BH q and the
                Bonferroni threshold computed WITHIN each garden x class), ->
                results/gemma_gwas/persite_gwas_{class}.npz + persite_gwas_summary.json

Usage:
    python gemma_persite_gwas.py scan --chrom chr1 --class snp
    python gemma_persite_gwas.py assemble
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
                              CHROMS, MIN_MAC_TEST, MIN_MAC_GRM, CALL_MIN)
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


def assemble():
    """Concatenate the per-chrom GEMMA parts and attach PER-GARDEN significance.

    No cross-site meta (scope decision 2026-08-25 -- the all-sites question belongs to the GEA
    track). Each garden is an independent scan, so multiple-testing correction is applied WITHIN
    a garden x class: BH q-values and the 0.05/M Bonferroni threshold. A marker significant in
    garden A and not garden B is NOT a tested contrast here -- comparing gardens is descriptive.
    Reported per garden: genomic-control lambda, Bonferroni count, BH count, best p.
    """
    order, Smat, sites, bio1 = _panel_order()
    S = len(sites)
    summary = {"n_founders": len(order), "n_gardens": S, "gardens": sites.tolist(),
               "bio1": bio1.tolist(), "scope": "per-garden only; no cross-site meta",
               "estimator": "GEMMA 0.98.5 -lmm 1 (Wald), univariate LMM, LOCO GRM via -k",
               "correction": "BH and Bonferroni computed within each garden x class"}
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
        # A marker is kept if it is testable in ANY garden (per-garden scans are independent);
        # untestable garden-cells stay NaN rather than dropping the whole marker.
        keep = np.isfinite(Zall).any(1)
        if (~keep).sum():
            print(f"[{clsname}] dropping {int((~keep).sum())} markers untestable in every garden",
                  flush=True)
        pos_all, chrom_all, Zall = pos_all[keep], chrom_all[keep], Zall[keep]
        M = len(pos_all)

        P = 2 * stats.norm.sf(np.abs(Zall))
        Q = np.full_like(P, np.nan)
        per_garden = []
        for si in range(S):
            ok = np.isfinite(P[:, si])
            m_ok = int(ok.sum()); bonf = 0.05 / max(m_ok, 1)
            Q[ok, si] = lib.bh(P[ok, si])
            per_garden.append(dict(
                garden=int(sites[si]), bio1=float(bio1[si]), n_markers=m_ok,
                lambda_gc=float(lib.lamgc(P[ok, si])),
                bonferroni_threshold=bonf,
                n_bonferroni=int((P[ok, si] < bonf).sum()),
                n_bh05=int((Q[ok, si] < 0.05).sum()),
                best_p=float(np.nanmin(P[ok, si])) if m_ok else float("nan")))
        np.savez(f"{OUT}/persite_gwas_{clsname}.npz", chrom=chrom_all, pos=pos_all,
                 Z=Zall, P=P, Q=Q, sites=sites, bio1=bio1,
                 lam_persite=np.array([g["lambda_gc"] for g in per_garden]))
        lam = np.array([g["lambda_gc"] for g in per_garden])
        nb = np.array([g["n_bonferroni"] for g in per_garden])
        summary[clsname] = dict(n_markers=M, per_garden=per_garden,
                                lambda_median=float(np.median(lam)),
                                lambda_min=float(lam.min()), lambda_max=float(lam.max()),
                                n_bonferroni_total=int(nb.sum()),
                                gardens_with_any_bonferroni=int((nb > 0).sum()))
        print(f"[{clsname}] M={M:,}  lambda median={np.median(lam):.3f} "
              f"({lam.min():.3f}-{lam.max():.3f})  Bonferroni hits={nb.sum()} "
              f"across {(nb>0).sum()}/{S} gardens", flush=True)

    json.dump(summary, open(f"{OUT}/persite_gwas_summary.json", "w"), indent=2, default=str)
    print(f"wrote {OUT}/persite_gwas_{{snp,nonsnp,sv}}.npz + persite_gwas_summary.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["scan", "assemble"])
    ap.add_argument("--chrom", choices=CHROMS)
    ap.add_argument("--class", dest="cls", choices=CLASSES)
    a = ap.parse_args()
    if a.stage == "scan":
        if not a.chrom or not a.cls:
            sys.exit("scan needs --chrom and --class")
        scan(a.chrom, a.cls)
    else:
        meta()
