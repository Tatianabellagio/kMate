#!/usr/bin/env python
"""Score every reviewed candidate on the evidence that made GPX6 stand out.

GPX6 is not strong because of one number. It passes a stack of independent checks, and
every one of them is now scripted. This runs the cheap ones for every reviewed candidate
(round-2 GEA own-axis + GWAS, round-1 survivors) so the one-by-one review starts from a
common table rather than from 340 figures:

  gradient      garden-level r (change from p0, averaged within each of the 31 gardens)
                against the candidate's own axis. Gardens, not pools, so the pool
                pseudo-replication that drives lambda is gone.
  pct_axis      where that |r| sits among 5,000 matched background indel/SV records ON
                THE SAME AXIS -- the number quoted for GPX6 (99.92nd).
  pct_sel       SELECTION-AWARE version: the candidate's axis was chosen as its best of
                22, so it is compared against each background record's best-of-22 |r|.
                This is the fair one; pct_axis flatters every candidate, GPX6 included.
  local_rank    rank of |r| among testable non-SNP records within +-25 kb (1 = lead)
  snp_rivals    SNPs within +-25 kb whose |r| is >= the candidate's (the "no kMate gain"
                count -- SNPs that mark the haplotype at least as well)
  region / tss  position relative to the target gene; signed distance upstream of TSS
  atac          bp of ATAC peak inside the REF footprint, peaks removed outright, tissues
  carriers      founder carriers / called, with SPLIT RECORDS MERGED (one SV event is
                often several same-size biallelic records; counting only the first
                undercounts -- the CRY2 bug)

Allele identity: every candidate enters through its reviewed store_row, which the
collider audit already resolved (memory `panel-multiallelic-pos-key-trap`). Nothing here
re-joins on chrom:pos.

Raw p throughout, no inflation correction (memory `raw-lfmm-over-wza-decision`).

env: kmate. Compute node, ~15-30 min. Writes results/candidate_evidence.csv.
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd, pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
import lib                                                       # noqa: E402
import axis_clusters as ac                                       # noqa: E402
import atac_overlap as AO                                        # noqa: E402

OUT = f"{HERE}/results"
PROJ = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
WIN, NBG, SEED = 25_000, 5000, 0


def candidates() -> pd.DataFrame:
    R = pd.read_csv(f"{OUT}/screen_visual_review_round2.csv")
    R = R[R.status == "done"]
    R = R[["set", "target_gene", "sym", "axis", "store_row", "verdict", "grade",
           "gwas_n_gardens"]]
    V = pd.read_csv(f"{OUT}/screen_visual_review.csv")
    V = V[~V.get("audit_2026_09_17", pd.Series("", index=V.index)).fillna("")
           .str.startswith("DROPPED")]
    # select BEFORE renaming: round-1 already carries its own `verdict` column, and
    # renaming final_tier onto it produced duplicate column labels
    V = V[["target_gene", "sym", "gea_best_axis", "store_row", "final_tier",
           "own_verdict"]].rename(columns={"gea_best_axis": "axis",
                                           "final_tier": "verdict", "own_verdict": "grade"})
    V["set"], V["gwas_n_gardens"] = "round1", np.nan
    V = V[["set", "target_gene", "sym", "axis", "store_row", "verdict", "grade",
           "gwas_n_gardens"]]
    C = pd.concat([R, V], ignore_index=True).drop_duplicates("store_row")
    return C.dropna(subset=["store_row"]).astype({"store_row": int}).reset_index(drop=True)


def site_means(M, p0, cols, site):
    X = np.asarray(M[:, cols], dtype=np.float32) - p0[cols][None, :]
    S = pd.DataFrame(X); S["site"] = site
    return S.groupby("site").mean()                    # sites x cols


def corr_axes(G, A):
    """|r| of every column of G (sites x vars) with every axis of A (sites x axes)."""
    A = A.reindex(G.index)
    Gv = G.values - np.nanmean(G.values, 0)
    out = {}
    for ax in A.columns:
        c = A[ax].values; ok = np.isfinite(c)
        cc = c[ok] - c[ok].mean()
        Z = Gv[ok]
        out[ax] = (np.nansum(Z * cc[:, None], 0)
                   / np.sqrt(np.nansum(Z ** 2, 0) * np.sum(cc ** 2)))
    return pd.DataFrame(out, index=G.columns)          # vars x axes (signed r)


def founder_support(vf, ch, pos, rl, al):
    """Carriers / called with same-size split records merged."""
    size = al - rl
    tol = max(10, int(abs(size) * 0.05)) if abs(size) > 50 else 0
    car, called = set(), set()
    for r in vf.fetch(ch, pos - 1, pos):
        if r.pos != pos:
            continue
        for a in r.alts or ():
            s = len(a) - len(r.ref)
            if np.sign(s) != np.sign(size) or abs(s - size) > tol:
                continue
            for smp in vf.header.samples:
                g = r.samples[smp].get("GT", (None,))[0]
                if g is not None:
                    called.add(smp)
                if g == 1:
                    car.add(smp)
    return len(car), len(called)


def main():
    rng = np.random.default_rng(SEED)
    C = candidates()
    print(f"{len(C)} candidates ({C.set.value_counts().to_dict()})", flush=True)

    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    K = pd.DataFrame({"chrom": idx["chrom"].astype(str), "pos": idx["pos"],
                      "ref_len": idx["ref_len"], "alt_len": idx["alt_len"]})
    C = C.join(K, on="store_row")
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    M = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy", mmap_mode="r")
    Ms = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_snp_af.npy", mmap_mode="r")
    p0s = np.load(f"{lib.AF_STORE}/p0_snp.npy")
    isn = np.load(f"{lib.AF_STORE}/index_snp.npz", allow_pickle=True)
    s_ch, s_pos = isn["chrom"].astype(str), isn["pos"]
    PM = pd.read_csv(ac.POOLMETA); site = PM.site.values
    A = ac.site_climate(); A.index = sorted(PM.site.astype(int).unique())

    # ---- background: same universe as the ATAC / audit backgrounds -----------------------
    B = pd.read_csv(f"{OUT}/repeat_variants_background_v2.csv.gz",
                    usecols=["chrom", "pos", "ref_len", "alt_len"])
    K["row"] = np.arange(len(K))
    brow = K.merge(B, on=["chrom", "pos", "ref_len", "alt_len"])["row"].values
    brow = np.sort(rng.choice(brow, NBG, replace=False))
    RB = corr_axes(site_means(M, p0, brow, site), A).abs()
    bg_axis = {ax: np.sort(RB[ax].dropna().values) for ax in A.columns}
    bg_sel = np.sort(RB.max(1).dropna().values)
    print(f"background: {len(brow):,} records; best-of-22 |r| median "
          f"{np.median(bg_sel):.3f}, p99 {np.percentile(bg_sel, 99):.3f}", flush=True)

    # ---- candidates: gradient on every axis ----------------------------------------------
    RC = corr_axes(site_means(M, p0, C.store_row.values, site), A)
    RC.index = C.index
    C["r_own"] = [RC.loc[i, ax] if ax in RC.columns else np.nan
                  for i, ax in zip(C.index, C.axis)]
    C["best_axis_all"] = RC.abs().idxmax(1)
    C["r_best_all"] = RC.abs().max(1)
    pct = lambda v, arr: 100 * np.searchsorted(arr, v) / len(arr)
    C["pct_axis"] = [pct(abs(r), bg_axis[ax]) if ax in bg_axis and np.isfinite(r) else np.nan
                     for r, ax in zip(C.r_own, C.axis)]
    C["pct_sel"] = [pct(v, bg_sel) for v in C.r_best_all]

    # ---- local rank + SNP rivals, on the own axis ----------------------------------------
    lr, sr, nloc, nsnp = [], [], [], []
    for c in C.itertuples():
        ax = c.axis if c.axis in A.columns else c.best_axis_all
        a = A[[ax]]
        sel = np.where((K.chrom.values == c.chrom) & (np.abs(K.pos.values - c.pos) <= WIN))[0]
        sel = sel[np.minimum(p0[sel], 1 - p0[sel]) > 0.05]
        sel = np.union1d(sel, [c.store_row])
        rn = corr_axes(site_means(M, p0, sel, site), a)[ax].abs()
        me = rn.iloc[int(np.where(sel == c.store_row)[0][0])]
        lr.append(int((rn > me + 1e-9).sum()) + 1); nloc.append(len(sel))
        ss = np.where((s_ch == c.chrom) & (np.abs(s_pos - c.pos) <= WIN))[0]
        ss = ss[np.minimum(p0s[ss], 1 - p0s[ss]) > 0.05]
        rs = corr_axes(site_means(Ms, p0s, ss, site), a)[ax].abs() if len(ss) else pd.Series([])
        sr.append(int((rs >= me - 1e-9).sum())); nsnp.append(len(ss))
    C["local_rank"], C["n_local"], C["snp_rivals"], C["n_snp_window"] = lr, nloc, sr, nsnp
    print("local ranks done", flush=True)

    # ---- position -------------------------------------------------------------------------
    VC = pd.read_csv(f"{OUT}/variants_classified.csv")
    VC["chrom"] = VC.chrom.str.replace("^chr", "Chr", regex=True)
    C = C.merge(VC[["chrom", "pos", "ref_len", "alt_len", "region", "tier"]]
                .drop_duplicates(["chrom", "pos", "ref_len", "alt_len"]),
                on=["chrom", "pos", "ref_len", "alt_len"], how="left")
    G = lib.load_genes().set_index("gene")
    tss_d = []
    for c in C.itertuples():
        if c.target_gene in G.index:
            g = G.loc[c.target_gene]
            g = g.iloc[0] if isinstance(g, pd.DataFrame) else g
            tss = g.end if g.strand == "-" else g.start
            tss_d.append((tss - c.pos) if g.strand == "-" else (c.pos - tss))
        else:
            tss_d.append(np.nan)
    C["dist_tss"] = tss_d        # negative = upstream of the TSS
    C["size"] = C.alt_len - C.ref_len

    # ---- ATAC -----------------------------------------------------------------------------
    peaks = AO.load_peaks()
    ov_bp, n_pk, whole, ntis = [], [], [], []
    for c in C.itertuples():
        ps, pe, num, _ = peaks[c.chrom]
        lo, hi = c.pos, c.pos + max(c.ref_len, 1) - 1
        h = np.where((pe >= lo) & (ps <= hi))[0]
        ov_bp.append(int(sum(min(pe[j], hi) - max(ps[j], lo) + 1 for j in h)))
        n_pk.append(len(h)); whole.append(int(sum((ps[j] >= lo) and (pe[j] <= hi) for j in h)))
        ntis.append(int(num[h].max()) if len(h) else 0)
    C["atac_bp"], C["atac_peaks"], C["atac_removed"], C["atac_tissues"] = ov_bp, n_pk, whole, ntis

    # ---- founder support ------------------------------------------------------------------
    car, cal = [], []
    vfs = {}
    for c in C.itertuples():
        f = f"{PROJ}/panel/arch3/{c.chrom.lower()}/merged_231_{c.chrom.lower()}_final.vcf.gz"
        vf = vfs.setdefault(f, pysam.VariantFile(f))
        a, b = founder_support(vf, c.chrom, int(c.pos), int(c.ref_len), int(c.alt_len))
        car.append(a); cal.append(b)
    C["carriers"], C["called"] = car, cal
    C["call_rate"] = (100 * C.called / 231).round(1)

    C.to_csv(f"{OUT}/candidate_evidence.csv", index=False)
    print(f"wrote {OUT}/candidate_evidence.csv ({len(C)} candidates)")


if __name__ == "__main__":
    main()
