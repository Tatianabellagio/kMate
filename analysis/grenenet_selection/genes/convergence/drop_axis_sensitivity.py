#!/usr/bin/env python
"""What survives if one climate axis is removed from the analysis entirely?

pc3 is the most inflated axis in the study (lambda 3.13), so candidates whose evidence
leans on it are structure-suspect. This asks the question directly instead of arguing it:
recompute everything that depends on the axis set, with that axis gone, and re-score.

Three things depend on the axis set, and each is recomputed here:

  1. POOL MEMBERSHIP. A variant is in the GEA pool because it is Bonferroni-significant on
     at least one axis. If its only significant axis is the dropped one and it is not a
     GWAS hit either, it would never have entered the pool -- it leaves.
  2. CLIMATE LINE. pct_sel is the best-of-22 garden-level |r| against a selection-aware
     background that is ALSO best-of-22; pct_quad is the same for the quadratic fit.
     Both sides lose the axis. Note the background gets easier too (a max over fewer axes
     is smaller), so a variant that never used pc3 can gain a little -- that is correct,
     not a bug, and is why the whole pool is recomputed rather than only the pc3 ones.
     Background rows are drawn exactly as evidence_matrix.py (seed 11) and
     nonlinear_gradient.py (seed 12) draw them, so the only difference is the axis.
  3. LOCUS LINE. local_rank and snp_rivals are computed on the variant's best axis. Where
     the best axis changes, they are recomputed on the new one (stage-b logic, +/-25 kb).

Lines that do not depend on the axis set -- GWAS, mechanism, chromatin, motif, expression --
are carried unchanged.

Usage: python drop_axis_sensitivity.py pc3          (one axis)
       python drop_axis_sensitivity.py pc3,bio15    (several at once)
Writes results/drop_<axes joined by +>_sensitivity.csv. env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", ".."))); sys.path.insert(0, HERE)
import lib, axis_clusters as ac                                  # noqa: E402
import candidate_evidence as CE                                  # noqa: E402
from nonlinear_gradient import quad                              # noqa: E402

OUT = f"{HERE}/results"
K4 = ["chrom", "pos", "ref_len", "alt_len"]


def background(K, seed):
    rng = np.random.default_rng(seed)
    Bg = pd.read_csv(f"{OUT}/repeat_variants_background_v2.csv.gz", usecols=K4)
    return np.sort(rng.choice(K.merge(Bg, on=K4)["row"].values, 4000, replace=False))


def local(v, axis, A, M, p0, kch, kpos, Ms, p0s, s_ch, s_pos, site):
    """stage-b local rank and SNP rivals, on a given axis."""
    a = A[[axis]]
    sel = np.where((kch == v.chrom) & (np.abs(kpos - v.pos) <= 25_000))[0]
    sel = np.union1d(sel[np.minimum(p0[sel], 1 - p0[sel]) > 0.05], [v.store_row])
    rn = CE.corr_axes(CE.site_means(M, p0, sel, site), a)[axis].abs()
    me = rn.iloc[int(np.where(sel == v.store_row)[0][0])]
    ss = np.where((s_ch == v.chrom) & (np.abs(s_pos - v.pos) <= 25_000))[0]
    ss = ss[np.minimum(p0s[ss], 1 - p0s[ss]) > 0.05]
    rs = (CE.corr_axes(CE.site_means(Ms, p0s, ss, site), a)[axis].abs()
          if len(ss) else pd.Series([], dtype=float))
    return int((rn > me + 1e-9).sum()) + 1, int((rs >= me - 1e-9).sum())


def main(drop_arg):
    drop = [a for a in drop_arg.split(",") if a]
    tag = "+".join(drop)
    V = pd.read_csv(f"{OUT}/evidence_matrix.csv")
    PM = pd.read_csv(ac.POOLMETA); site = PM.site.values
    A = ac.site_climate(); A.index = sorted(PM.site.astype(int).unique())
    bad = [a for a in drop if a not in A.columns]
    assert not bad, f"{bad} not axes: {list(A.columns)}"
    A2 = A.drop(columns=drop)

    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    K = pd.DataFrame({"chrom": idx["chrom"].astype(str), "pos": idx["pos"],
                      "ref_len": idx["ref_len"], "alt_len": idx["alt_len"]})
    K["row"] = np.arange(len(K))
    M = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy", mmap_mode="r")
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")

    # 1. pool membership
    sig = V.gea_sig_axes.fillna("").astype(str).str.split(",")
    others = sig.map(lambda s: [a for a in s if a and a not in drop])
    V["pool_kept"] = others.map(len).gt(0) | V.in_gwas.astype(bool)

    # 2. climate line, both sides without the axis
    Gb = CE.site_means(M, p0, background(K, 11), site)
    bg_sel = np.sort(CE.corr_axes(Gb, A2).abs().max(axis=1).values)
    Gq = CE.site_means(M, p0, background(K, 12), site)
    _, bg_r2, _, _ = quad(Gq, A2); bg_r2 = np.sort(bg_r2)
    Gc = CE.site_means(M, p0, V.store_row.values, site); Gc.columns = V.index
    RC = CE.corr_axes(Gc, A2)
    V["best_axis_2"] = RC.abs().idxmax(axis=1).values
    V["pct_sel_2"] = 100 * np.searchsorted(bg_sel, RC.abs().max(axis=1).values) / len(bg_sel)
    _, r2, _, _ = quad(Gc, A2)
    V["pct_quad_2"] = 100 * np.searchsorted(bg_r2, r2) / len(bg_r2)
    V["L_climate_2"] = ((V.pct_sel_2 >= 97) | (V.pct_quad_2 >= 97)).astype(int)

    # 3. locus line, recomputed only where the best axis moved
    V["local_rank_2"], V["snp_rivals_2"] = V.local_rank, V.snp_rivals
    # also for the leavers that were candidates (n_lines >= 3): "would it still be a good
    # candidate on the honest axes?" is exactly the question for them, even though the
    # GEA scan alone would no longer have nominated them
    moved = V.index[(V.best_axis_2 != V.best_axis) & (V.pool_kept | (V.n_lines >= 3))]
    print(f"best axis changes for {len(moved)} variants (pool-kept or candidates); "
          f"recomputing local rank", flush=True)
    if len(moved):
        kch, kpos = idx["chrom"].astype(str), idx["pos"]
        Ms = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_snp_af.npy", mmap_mode="r")
        p0s = np.load(f"{lib.AF_STORE}/p0_snp.npy")
        isn = np.load(f"{lib.AF_STORE}/index_snp.npz", allow_pickle=True)
        s_ch, s_pos = isn["chrom"].astype(str), isn["pos"]
        for n, i in enumerate(moved, 1):
            lr, sr = local(V.loc[i], V.loc[i, "best_axis_2"], A2, M, p0, kch, kpos,
                           Ms, p0s, s_ch, s_pos, site)
            V.loc[i, "local_rank_2"], V.loc[i, "snp_rivals_2"] = lr, sr
            if n % 50 == 0:
                print(f"  {n}/{len(moved)}", flush=True)
    V["L_locus_2"] = (V.local_rank_2 <= 2).astype(int)

    kept_lines = ["L_gwas", "L_mechanism", "L_chromatin", "L_motif", "L_expression"]
    # n_lines_2 is scored for everyone; pool_kept says whether the GEA scan would still have
    # nominated the variant. Read the two together: a leaver with a high n_lines_2 is a gene
    # the honest axes support but the scan would not have found on its own.
    V["n_lines_2"] = V[["L_climate_2", "L_locus_2"] + kept_lines].sum(axis=1)
    V.to_csv(f"{OUT}/drop_{tag}_sensitivity.csv", index=False)

    print(f"\ndropping {tag}: pool {len(V)} -> {int(V.pool_kept.sum())} "
          f"({int((~V.pool_kept).sum())} leave: significant only on {tag}, not GWAS hits)")
    print("n_lines before:", V.n_lines.value_counts().sort_index().to_dict())
    print("n_lines after (pool-kept):", V[V.pool_kept].n_lines_2.value_counts().sort_index().to_dict())


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "pc3")
