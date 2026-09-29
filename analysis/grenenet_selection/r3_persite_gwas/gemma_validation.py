#!/usr/bin/env python
"""Validate the in-house per-site LOCO-EMMAX scan (class_split_gwas.py) against GEMMA.

Motivation (user, 2026-08-25): the per-site association model is currently our own ~70-line
EMMAX/P3D implementation (founder_genotype.emma_reml_delta + a closed-form 2x2 GLS in
class_split_gwas.py). Reviewers trust a named, published tool more than a hand-rolled one, so
we either swap in GEMMA or show ours reproduces it. This script measures the agreement.

What it does, on a SUBSET (one chrom x one class x a few sites, subsampled markers):
  1. Rebuilds the EXACT pooled per-chrom GRM pass-1 of class_split_gwas.py (MAC>=12, call>=90%,
     all classes pooled) by importing its own helpers -- so the kinship, marker filters, founder
     order and missing-imputation are identical by construction, not re-derived.
  2. Forms the LOCO GRM for the pilot chromosome: (sum_chrom K - K_chrom) / (Mtot - m_chrom).
  3. Writes GEMMA inputs: BIMBAM mean-genotype (dosages as given, 0/1 haploid coding, missing
     already imputed to the marker mean exactly as our scan does), the rank-normalized per-site
     phenotype matrix, and the LOCO GRM as GEMMA's -k relatedness matrix.
  4. Runs `gemma -lmm 1` (Wald) per site with GEMMA's own QC filters DISABLED (-maf/-miss/-r2/-hwe
     wide open) so GEMMA tests exactly our marker set and no marker-count difference confounds
     the comparison.
  5. Compares GEMMA's z = beta/se against our stored per-site Z (class_gwas_<class>.npz, matched
     by position): Pearson/Spearman on z, on -log10 p, and the max absolute discrepancy.

The expected disagreement is small but NOT zero, and it is one-directional by design: GEMMA's
-lmm 1 re-estimates the variance component per marker (exact Wald), while ours fixes delta under
the null per (site, chrom) -- the standard EMMAX/P3D approximation. So GEMMA is the stricter
model and ours should track it closely, diverging most at the largest-effect markers.

Writes results/gemma_validation/{pooled_grms.npz, gemma_io/, concordance.csv, summary.json}.
Env: needs BOTH `kmate` (python side) and the `gemma` binary (env gwas_tools). Run as:
    conda activate kmate && PATH=$HOME/miniforge3/envs/gwas_tools/bin:$PATH python gemma_validation.py
Compute-node only (dense per-chrom marker blocks; pass 1 peaks ~2-3 GB).
"""
from __future__ import annotations
import os, sys, json, time, shutil, subprocess
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
# reuse the production helpers verbatim so filters/order/imputation cannot drift
from class_split_gwas import (_panel_order, _load_chrom_212, _classes, _dense_imputed, _ZZt,
                              CHROMS, MIN_MAC_TEST, MIN_MAC_GRM, CALL_MIN)

VAREXP = f"{lib.GEA}/r3_persite_gwas/results/varexp"
OUT = f"{lib.GEA}/r3_persite_gwas/results/gemma_validation"
IO = f"{OUT}/gemma_io"

PILOT_CHROM = "chr1"
PILOT_CLASS = "snp"          # snp | nonsnp | sv
N_MARKERS = 20_000           # subsample for the pilot; enough to pin a correlation
N_SITES = 3
SEED = 0

# GEMMA 0.98.5 BUG (measured here 2026-08-25, not a property of our data): any single run with
# >= 20,000 markers dies with "Enforce failed for N>0 in src/fastblas.cpp ... fast_cblas_dgemm"
# AFTER completing the scan but BEFORE writing .assoc.txt. Verified boundary: 19,999 markers OK,
# 20,000 / 20,001 / 21,000 / 30,000 / 45,000 all fail; independent of OPENBLAS_NUM_THREADS and of
# the QC flags. GEMMA's fast-BLAS path handles exactly one 20k block and mishandles the next.
# Workaround: split the marker set into chunks below the limit and concatenate the per-chunk
# association output. Cheap -- the per-run cost GEMMA repeats is a 231x231 eigendecomposition.
GEMMA_MAX = 19_999


def build_pooled_grms(order, N):
    """Pass 1 of class_split_gwas.py: per-chrom pooled GRM + genome total. Cached."""
    cache = f"{OUT}/pooled_grms.npz"
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return {c: (z[f"K_{c}"], int(z[f"m_{c}"])) for c in CHROMS}, z["Ksum"], int(z["Mtot"])
    t0 = time.time()
    Kchrom, Ksum, Mtot = {}, np.zeros((N, N)), 0
    for cl in CHROMS:
        pos, rl, al, vp, vc, n_alt, n_cal = _load_chrom_212(cl, order)
        mac = np.minimum(n_alt, N - n_alt)
        gcols = np.where((mac >= MIN_MAC_GRM) & (n_cal >= CALL_MIN * N))[0]
        g, p = _dense_imputed(vp, vc, n_alt, n_cal, gcols)
        K = _ZZt(g, p)
        del g
        Kchrom[cl] = (K, len(gcols)); Ksum += K; Mtot += len(gcols)
        print(f"  grm {cl}: {len(gcols):,} markers ({time.time()-t0:.0f}s)", flush=True)
    np.savez(cache, Ksum=Ksum, Mtot=Mtot,
             **{f"K_{c}": Kchrom[c][0] for c in CHROMS},
             **{f"m_{c}": Kchrom[c][1] for c in CHROMS})
    return Kchrom, Ksum, Mtot


def main():
    os.makedirs(IO, exist_ok=True)
    t0 = time.time()
    if shutil.which("gemma") is None:
        sys.exit("gemma not on PATH -- prepend $HOME/miniforge3/envs/gwas_tools/bin")

    order, Smat, sites, bio1 = _panel_order()
    N = len(order)
    print(f"{N} founders, {len(sites)} sites; pilot = {PILOT_CHROM}/{PILOT_CLASS}, "
          f"{N_MARKERS:,} markers x {N_SITES} sites", flush=True)

    Kchrom, Ksum, Mtot = build_pooled_grms(order, N)
    K_loco = (Ksum - Kchrom[PILOT_CHROM][0]) / max(Mtot - Kchrom[PILOT_CHROM][1], 1)

    # ---- pilot test markers: identical filter to the production scan ----
    pos, rl, al, vp, vc, n_alt, n_cal = _load_chrom_212(PILOT_CHROM, order)
    mac = np.minimum(n_alt, N - n_alt)
    vclass = _classes(rl, al)
    keep = (mac >= MIN_MAC_TEST) & (n_cal >= CALL_MIN * N)
    want = {"snp": vclass == 0, "nonsnp": vclass != 0, "sv": vclass == 2}[PILOT_CLASS]
    idx_all = np.where(keep & want)[0]
    rng = np.random.default_rng(SEED)
    idx = np.sort(rng.choice(idx_all, size=min(N_MARKERS, len(idx_all)), replace=False))
    g, p = _dense_imputed(vp, vc, n_alt, n_cal, idx)      # [N founders x M markers]
    print(f"  {len(idx_all):,} eligible {PILOT_CLASS} markers on {PILOT_CHROM}, "
          f"using {len(idx):,} ({time.time()-t0:.0f}s)", flush=True)

    # ---- GEMMA inputs ----
    # BIMBAM mean genotype: rs, minor, major, then one dosage per founder
    marker_pos = pos[idx]
    # rs carries the marker's ORDINAL in the pilot set, not just chrom:pos. Keying the
    # comparison on chrom:pos is WRONG here: 2.14% of chr1 positions carry >1 biallelic record
    # (multiallelic decomposition in arch3), so a position key silently matches a different alt
    # allele and manufactures sign flips. Measured 2026-08-25: it dragged pearson(z) from 0.9975
    # (correct matching) down to 0.9848, with 29% sign disagreement on the affected 2.2%.
    rs = np.array([f"{PILOT_CHROM}:{q}:{j}" for j, q in enumerate(marker_pos)])
    gt = g.T                                              # [M x N]
    chunks = np.array_split(np.arange(len(idx)), int(np.ceil(len(idx) / GEMMA_MAX)))
    for ci, ch in enumerate(chunks):                      # see GEMMA_MAX: >=20k per run crashes
        with open(f"{IO}/geno_{ci}.bimbam", "w") as fh:
            for j in ch:
                fh.write(f"{rs[j]},A,T," + ",".join(f"{v:.6g}" for v in gt[j]) + "\n")
    print(f"  geno split into {len(chunks)} chunk(s) of <= {max(len(c) for c in chunks):,}",
          flush=True)
    site_use = list(range(min(N_SITES, len(sites))))
    P = np.column_stack([stats.norm.ppf((stats.rankdata(Smat[si]) - 0.5) / N) for si in site_use])
    np.savetxt(f"{IO}/pheno.txt", P, fmt="%.10f")
    np.savetxt(f"{IO}/kin.txt", K_loco, fmt="%.10f")
    print(f"  wrote GEMMA inputs ({time.time()-t0:.0f}s)", flush=True)

    # ---- run GEMMA per site (QC filters wide open: test exactly our marker set) ----
    ours = np.load(f"{VAREXP}/class_gwas_{PILOT_CLASS}.npz", allow_pickle=True)
    ours_mask = ours["chrom"].astype(str) == PILOT_CHROM
    ours_pos = ours["pos"][ours_mask]
    ours_Z = ours["Z"][ours_mask]
    # Map each pilot marker to its production row by ORDER, not position. The stored npz is
    # pos[idx_all] minus the handful of markers the near-singular guard dropped, so ours_pos is
    # an in-order subsequence of pos_full -- recover the alignment with a two-pointer walk.
    pos_full = pos[idx_all]
    row_of_full = np.full(len(pos_full), -1, np.int64)
    j = 0
    for i_, q in enumerate(pos_full):
        if j < len(ours_pos) and ours_pos[j] == q:
            row_of_full[i_] = j; j += 1
    assert j == len(ours_pos), f"row alignment failed ({j} matched of {len(ours_pos)})"
    n_drop = int((row_of_full < 0).sum())
    row_for_pilot = row_of_full[np.searchsorted(idx_all, idx)]
    print(f"  aligned production rows; {n_drop} marker(s) dropped by the near-singular guard",
          flush=True)

    rows, summary = [], {"chrom": PILOT_CHROM, "class": PILOT_CLASS,
                         "n_markers_requested": int(len(idx)), "sites": []}
    for k, si in enumerate(site_use):
        parts = []
        for ci in range(len(chunks)):
            pref = f"site{sites[si]}_c{ci}"
            cmd = ["gemma", "-g", f"{IO}/geno_{ci}.bimbam", "-p", f"{IO}/pheno.txt",
                   "-k", f"{IO}/kin.txt", "-lmm", "1", "-n", str(k + 1),
                   "-maf", "0", "-miss", "1", "-r2", "1", "-hwe", "0", "-outdir", IO, "-o", pref]
            r = subprocess.run(cmd, capture_output=True, text=True)
            assoc = f"{IO}/{pref}.assoc.txt"
            if not os.path.exists(assoc):
                print(r.stdout[-2000:]); print(r.stderr[-2000:])
                sys.exit(f"GEMMA failed for {pref}")
            parts.append(np.genfromtxt(assoc, names=True, dtype=None, encoding=None))
        a = np.concatenate(parts) if len(parts) > 1 else parts[0]
        gj = np.array([int(str(s).split(":")[2]) for s in a["rs"]])   # pilot ordinal, exact
        gz = a["beta"] / a["se"]
        rr = row_for_pilot[gj]
        keep_m = rr >= 0
        gpos = marker_pos[gj]
        oz = ours_Z[rr[keep_m], si]
        gz_m = gz[keep_m]
        ok = np.isfinite(oz) & np.isfinite(gz_m)
        pear = float(np.corrcoef(oz[ok], gz_m[ok])[0, 1])
        spear = float(stats.spearmanr(oz[ok], gz_m[ok]).statistic)
        onlp = -np.log10(np.clip(2 * stats.norm.sf(np.abs(oz[ok])), 1e-300, 1))
        gnlp = -np.log10(np.clip(a["p_wald"][keep_m][ok], 1e-300, 1))
        summary["sites"].append(dict(
            site=int(sites[si]), n_matched=int(ok.sum()),
            pearson_z=pear, spearman_z=spear,
            pearson_nlogp=float(np.corrcoef(onlp, gnlp)[0, 1]),
            max_abs_z_diff=float(np.abs(oz[ok] - gz_m[ok]).max()),
            median_abs_z_diff=float(np.median(np.abs(oz[ok] - gz_m[ok]))),
            our_lambda=float(np.median(oz[ok] ** 2) / stats.chi2.ppf(0.5, 1)),
            gemma_lambda=float(np.median(gz_m[ok] ** 2) / stats.chi2.ppf(0.5, 1))))
        for q, a_, b_ in zip(gpos[keep_m][ok], oz[ok], gz_m[ok]):
            rows.append((int(sites[si]), int(q), float(a_), float(b_)))
        s = summary["sites"][-1]
        print(f"  site {sites[si]}: n={s['n_matched']:,} pearson(z)={pear:.6f} "
              f"spearman={spear:.6f} max|dz|={s['max_abs_z_diff']:.3g} "
              f"lambda ours={s['our_lambda']:.3f} gemma={s['gemma_lambda']:.3f} "
              f"({time.time()-t0:.0f}s)", flush=True)

    import csv
    with open(f"{OUT}/concordance.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["site", "pos", "z_ours", "z_gemma"]); w.writerows(rows)
    allp = [s["pearson_z"] for s in summary["sites"]]
    summary["min_pearson_z"] = float(min(allp))
    summary["mean_pearson_z"] = float(np.mean(allp))
    json.dump(summary, open(f"{OUT}/summary.json", "w"), indent=2)
    print(f"\nmin pearson(z) across sites = {min(allp):.6f}")
    print(f"wrote {OUT}/concordance.csv + summary.json ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
