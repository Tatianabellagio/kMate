#!/usr/bin/env python
"""Why is lambda so much higher on some climate axes, and is it a low-frequency effect?

Two questions the axis-by-axis lambda table (`raw_block_significance_lfmm.csv`) raises but
cannot answer, asked here directly. Nothing is corrected -- lambda is measured, never
divided out (memory `raw-lfmm-over-wza-decision`).

  --part maf     Is the inflation carried by variants near the MAF 5% floor?
                 Per (class x axis): lambda and Bonferroni hit rate WITHIN MAF strata, and
                 a floor sweep (0.05 / 0.10 / 0.15 / 0.20 / 0.25) recomputing lambda, the
                 threshold and the hit count on each retained set. This is the test that
                 decides whether raising --maf-min would buy calibration.
                 -> results/multiaxis/axis_inflation_maf.csv

  --part struct  Does an axis's lambda track how much of the axis lines up with the
                 dominant structure of the response? The LFMM response is Y = dp
                 [352 pools x records] and the model absorbs K = 16 latent factors in POOL
                 space. So decompose pool space itself: eigenvectors of the pool x pool
                 Gram matrix of the (column-standardised) dp, i.e. the PCs of the
                 genome-wide allele-frequency-change pattern, then ask what share of each
                 climate axis's variance sits on PCs 1-16 (the part the model can absorb)
                 versus further out (the part it cannot).
                 -> results/multiaxis/axis_inflation_struct.csv, ..._poolpcs.csv

  --part gardens Same residual covariation, but leaving out the four founder-sweep gardens
                 (4, 43, 45, 32 -- the set the candidate review already flags as a WEAK
                 pattern). bio15's alignment sits on the pool PC those gardens load, so this
                 asks whether its inflation IS those gardens. Everything is recomputed on
                 the retained pools, including the latent-factor basis.
                 -> results/multiaxis/axis_inflation_gardens.csv

Note on the unit: all 31 gardens were sown from the SAME seed mix, so there is no ancestral
population structure across gardens in the usual GEA sense. The "structure" the latent
factors face is the shared garden-level RESPONSE -- founder lineages rising and falling
together -- plus the fact that ~11 pools share each garden's single climate value (352 rows,
31 distinct X values; METHODS_lfmm.md §1, §7). The pseudo-replication is identical on every
axis, so it cannot explain the ranking; the alignment measured here can.

env: kmate. Compute node (--part struct reads the 3.2 GB gen-9 non-SNP AF matrix once).
"""
from __future__ import annotations
import argparse
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..")))
import lib                                                       # noqa: E402

MULTI = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis"
WZAIN = f"{MULTI}/tile_cr50"   # production (events, MAF 0.05, call rate 0.50); was wza_in_clq09_tile until 2026-10-04
LAMTAB = f"{MULTI}/raw_block_significance_lfmm_tile_cr50.csv"
POOLMETA = f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp.meta.csv"
AFMAT = f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy"
P0 = f"{lib.GEA}/common/results/af_store/p0_nonsnp.npy"

BIOS = [f"bio{i}" for i in range(1, 20)]
AXES = BIOS + ["pc1", "pc2", "pc3"]
CLASSES = ["snp", "nonsnp", "smallindel", "sv"]
MAF_MIN = 0.05                      # the production floor
MAF_BINS = [0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.50001]
FLOORS = [0.05, 0.10, 0.15, 0.20, 0.25]
K_LFMM = 16                         # latent factors the production runs absorb
SWEEP_GARDENS = [4, 43, 45, 32]     # the founder-sweep set (`r_garden4_set` in the review)
N_SUB = 400_000                     # records sampled for the pool-space work (of 2.25 M)


def lam_of(p: np.ndarray) -> float:
    """Genomic inflation factor from p-values: median chi2(1 df) / 0.4549364."""
    return float(np.median(chi2.isf(np.clip(p, 1e-300, 1.0), 1)) / 0.4549364)


# ---------------------------------------------------------------------------- part: maf
def part_maf(classes: list[str]) -> None:
    rows_s, rows_f = [], []
    for cls in classes:
        for ax in AXES:
            f = f"{WZAIN}/lfmm_{cls}_gen9_{ax}.csv"
            if not os.path.exists(f):
                print(f"  MISSING {os.path.basename(f)}", flush=True)
                continue
            d = pd.read_csv(f, usecols=["MAF", "pval"])
            d = d[d.MAF > MAF_MIN]
            thr0 = 0.05 / len(d)                       # the production threshold
            # --- lambda and hit rate within MAF strata, at the production threshold
            b = pd.cut(d.MAF, MAF_BINS, right=False)
            for k, g in d.groupby(b, observed=True):
                rows_s.append(dict(cls=cls, axis=ax, maf_lo=k.left, maf_hi=k.right,
                                   n=len(g), lam=lam_of(g.pval.values),
                                   hits=int((g.pval < thr0).sum()),
                                   hits_per_10k=1e4 * float((g.pval < thr0).mean())))
            # --- floor sweep: everything recomputed on the retained set
            sig0 = d.pval < thr0
            for fl in FLOORS:
                k = d.MAF >= fl if fl > MAF_MIN else d.MAF > MAF_MIN
                sub = d[k]
                thr = 0.05 / len(sub)
                rows_f.append(dict(
                    cls=cls, axis=ax, floor=fl, n=len(sub),
                    frac_records_kept=len(sub) / len(d),
                    lam=lam_of(sub.pval.values), bonf=thr,
                    hits=int((sub.pval < thr).sum()),
                    hits_old_thr=int((sub.pval < thr0).sum()),
                    hits_kept_of_orig=int((sig0 & k).sum()),
                    hits_orig=int(sig0.sum())))
            print(f"  {cls:11s} {ax:5s} n={len(d):>9,} lam={lam_of(d.pval.values):.2f} "
                  f"hits={int(sig0.sum()):>5,}", flush=True)
        # written per class: the snp scan alone is ~45 min, so a run that dies partway
        # keeps every class it finished
        pd.DataFrame([r for r in rows_s if r["cls"] == cls]).to_csv(
            f"{MULTI}/axis_inflation_maf_{cls}.csv", index=False)
        pd.DataFrame([r for r in rows_f if r["cls"] == cls]).to_csv(
            f"{MULTI}/axis_inflation_floor_{cls}.csv", index=False)
        print(f"  -> wrote axis_inflation_{{maf,floor}}_{cls}.csv", flush=True)
    for tag, rows in (("maf", rows_s), ("floor", rows_f)):
        have = [f"{MULTI}/axis_inflation_{tag}_{c}.csv" for c in CLASSES]
        if all(os.path.exists(f) for f in have):     # all four classes on disk: combine
            pd.concat([pd.read_csv(f) for f in have]).to_csv(
                f"{MULTI}/axis_inflation_{tag}.csv", index=False)
            print(f"wrote {MULTI}/axis_inflation_{tag}.csv (all classes)")


# ------------------------------------------------------------------------- part: struct
def pool_axes() -> pd.DataFrame:
    """352 pools x 22 axes, built exactly as build_class_matrices.py builds them."""
    d = pd.read_csv(POOLMETA)
    m = d[BIOS].to_numpy(float)
    Z = (m - m.mean(0)) / m.std(0, ddof=0)
    U, S, Vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
    for k in range(3):
        sc = U[:, k] * S[k]
        flip = (np.corrcoef(sc, m[:, 0])[0, 1] < 0 if k == 0
                else Vt[k][np.argmax(np.abs(Vt[k]))] < 0)
        if flip:
            sc = -sc
        d[f"pc{k + 1}"] = (sc - sc.mean()) / sc.std(ddof=0)
    return d


def pool_gram(chunk: int = 40_000) -> tuple[np.ndarray, int]:
    """Gram matrix of the column-standardised dp over the records LFMM tested.

    Column-standardised (each record's dp scaled to unit SD across pools) because the
    per-record test statistic is scale-free: what matters is the DIRECTION in pool space
    each record's response points, not its amplitude. MAF is recomputed here the way
    build_class_matrices.py defines it (mean gen-9 AF across pools, folded).
    """
    af = np.load(AFMAT, mmap_mode="r")
    p0 = np.load(P0)
    n_pools, n_rec = af.shape
    G = np.zeros((n_pools, n_pools))
    kept = 0
    for i in range(0, n_rec, chunk):
        a = np.asarray(af[:, i:i + chunk], dtype=np.float64)
        pbar = a.mean(0)
        maf = np.minimum(pbar, 1.0 - pbar)
        k = maf >= MAF_MIN
        if not k.any():
            continue
        y = a[:, k] - p0[i:i + chunk][k]               # dp, as build_lfmm_input.py writes it
        y -= y.mean(0)                                 # LFMM centres Y
        sd = y.std(0, ddof=0)
        good = sd > 0
        y = y[:, good] / sd[good]
        G += y @ y.T
        kept += int(good.sum())
        if (i // chunk) % 10 == 0:
            print(f"  {i:>10,}/{n_rec:,}  kept {kept:,}", flush=True)
    return G / kept, kept


def part_struct() -> None:
    A = pool_axes()
    G, kept = pool_gram()
    print(f"\nGram over {kept:,} records (LFMM tested 703,789 non-SNP)")
    w, V = np.linalg.eigh(G)
    w, V = w[::-1], V[:, ::-1]                         # descending
    vexp = w / w.sum()

    rows = []
    for ax in AXES:
        x = A[ax].to_numpy(float)
        x = x - x.mean()
        x = x / np.linalg.norm(x)
        r2 = (V.T @ x) ** 2                            # share of the axis on each pool PC
        best = int(np.argmax(r2))
        rows.append(dict(
            axis=ax,
            r2_pc1=r2[0], r2_pc1_5=r2[:5].sum(), r2_in_K=r2[:K_LFMM].sum(),
            r2_beyond_K=r2[K_LFMM:].sum(),
            best_pc=best + 1, r2_best_pc=r2[best],
            # where the axis's alignment sits on the PC ranking: 1 = all on PC1
            pc_rank_centroid=float((np.arange(1, len(r2) + 1) * r2).sum() / r2.sum()),
            var_at_best_pc=vexp[best]))
    S = pd.DataFrame(rows)

    lam = pd.read_csv(LAMTAB)
    lam = lam[(lam.model == "lfmm") & (lam.cls == "nonsnp")][["axis", "lam"]]
    S = S.merge(lam, on="axis", how="left").sort_values("lam", ascending=False)
    S.to_csv(f"{MULTI}/axis_inflation_struct.csv", index=False)
    pd.DataFrame(dict(pc=np.arange(1, len(w) + 1), eigenvalue=w,
                      var_explained=vexp,
                      cum_var=np.cumsum(vexp))).to_csv(
        f"{MULTI}/axis_inflation_poolpcs.csv", index=False)

    pd.set_option("display.width", 200)
    print("\nvariance of the dp pool-space structure: "
          f"PC1 {vexp[0]:.1%}, PC1-5 {vexp[:5].sum():.1%}, "
          f"PC1-{K_LFMM} {vexp[:K_LFMM].sum():.1%}")
    print(S.round(4).to_string(index=False))
    print(f"\nSpearman lam vs r2_in_K: {S.lam.corr(S.r2_in_K, method='spearman'):.3f}"
          f"   vs r2_pc1: {S.lam.corr(S.r2_pc1, method='spearman'):.3f}"
          f"   vs pc_rank_centroid: "
          f"{S.lam.corr(S.pc_rank_centroid, method='spearman'):.3f}")
    print(f"\nwrote {MULTI}/axis_inflation_struct.csv, axis_inflation_poolpcs.csv")


# ------------------------------------------------------------------------ part: gardens
def dp_subsample(seed: int = 7, chunk: int = 40_000) -> np.ndarray:
    """dp for a random N_SUB-record sample of what LFMM tested, [352 pools x records].

    A sample, not the whole matrix: the quantities below are genome-wide MEANS over
    records, so ~125k tested records pin them to three decimals.

    Read in CONTIGUOUS column chunks, selecting the sampled columns within each chunk. One
    fancy index of 400k scattered columns over the 3.2 GB memmap does not finish in 20
    minutes on shared scratch (it is millions of strided reads -- the hazard
    build_class_matrices.py:88 warns about); chunked, the same data takes about a minute.
    """
    af = np.load(AFMAT, mmap_mode="r")
    p0 = np.load(P0)
    n_rec = af.shape[1]
    rng = np.random.default_rng(seed)
    cols = np.sort(rng.choice(n_rec, N_SUB, replace=False))
    out = []
    for i in range(0, n_rec, chunk):
        c = cols[(cols >= i) & (cols < i + chunk)]
        if not c.size:
            continue
        a = np.asarray(af[:, i:i + chunk], dtype=np.float64)[:, c - i]
        pbar = a.mean(0)
        keep = np.minimum(pbar, 1.0 - pbar) >= MAF_MIN
        if keep.any():
            out.append(a[:, keep] - p0[c][keep])
    return np.concatenate(out, axis=1)


def resid_mean_r2(Y0: np.ndarray, A: pd.DataFrame, mask: np.ndarray,
                  K: int = K_LFMM) -> tuple[dict[str, float], float]:
    """Mean r2(dp, axis) over records AFTER projecting out the top-K pool-space PCs.

    The latent-factor basis is recomputed on the retained pools, so a leave-out is a
    genuine refit of the correction, not the full-data basis reapplied. Returns the
    per-axis means and the noise floor 1/(n-1-K) they should be read against.
    """
    Y = Y0[mask]
    Y = Y - Y.mean(0)
    sd = Y.std(0, ddof=0)
    g = sd > 0
    ys = Y[:, g] / sd[g]
    G = (ys @ ys.T) / int(g.sum())
    _, V = np.linalg.eigh(G)
    P = V[:, ::-1][:, :K]
    Yr = ys - P @ (P.T @ ys)
    out = {}
    for ax in AXES:
        x = A[ax].to_numpy(float)[mask]
        xc = x - x.mean()
        xr = xc - P @ (P.T @ xc)
        num = Yr.T @ xr
        den = np.sqrt((Yr ** 2).sum(0) * (xr @ xr))
        r = np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0)
        out[ax] = float(np.mean(r ** 2))
    n = int(mask.sum())
    return out, 1.0 / (n - 1 - K)


def part_gardens() -> None:
    A = pool_axes()
    Y0 = dp_subsample()
    print(f"dp subsample: {Y0.shape[1]:,} tested records x {Y0.shape[0]} pools", flush=True)
    full = np.ones(len(A), bool)
    drop = A.site.isin(SWEEP_GARDENS).to_numpy()
    r_all, f_all = resid_mean_r2(Y0, A, full)
    print(f"  full: noise floor {f_all:.5f}", flush=True)
    r_out, f_out = resid_mean_r2(Y0, A, ~drop)
    print(f"  without gardens {SWEEP_GARDENS}: noise floor {f_out:.5f}", flush=True)

    lam = pd.read_csv(LAMTAB)
    lam = lam[(lam.model == "lfmm") & (lam.cls == "nonsnp")].set_index("axis").lam
    T = pd.DataFrame({
        "lam": lam,
        "resid_r2": pd.Series(r_all), "x_floor": pd.Series(r_all) / f_all,
        "resid_r2_no_sweep": pd.Series(r_out), "x_floor_no_sweep": pd.Series(r_out) / f_out,
    })
    T["change"] = T.x_floor_no_sweep / T.x_floor
    T.index.name = "axis"
    T = T.sort_values("lam", ascending=False)
    T.to_csv(f"{MULTI}/axis_inflation_gardens.csv")
    pd.set_option("display.width", 200)
    print(f"\nn pools {int(full.sum())} -> {int((~drop).sum())} "
          f"({int(drop.sum())} pools in gardens {SWEEP_GARDENS})")
    print("x_floor = residual covariation in units of the noise floor; 1.0 <=> lambda 1")
    print(T.round(3).to_string())
    print(f"\nwrote {MULTI}/axis_inflation_gardens.csv")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["maf", "struct", "gardens"], required=True)
    ap.add_argument("--classes", default=",".join(CLASSES))
    a = ap.parse_args()
    if a.part == "maf":
        part_maf(a.classes.split(","))
    elif a.part == "struct":
        part_struct()
    else:
        part_gardens()
