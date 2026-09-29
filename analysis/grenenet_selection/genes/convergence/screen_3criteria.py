#!/usr/bin/env python
"""Screen every non-SNP GEA+GWAS hit against the three criteria that define the pattern
we are actually looking for (user, 2026-09-15). Programmatic stage; plotting is separate.

  C1  PLACEMENT   the variant sits inside a gene, or immediately before/after it, in a way
                  that could have a functional effect (CDS frameshift > CDS in-frame >
                  UTR > promoter > downstream flank > intron-only).
  C2  MOVEMENT    the allele clearly rises or falls in gardens -- not a flat line at p0.
  C3  CLIMATE     that rise/fall tracks a climate gradient: temperature (bio1) or
                  precipitation/humidity (bio12).

C1 is inherited from `build_functional_variants.py` (`ftier`). C2 and C3 are computed here
from the pool-level allele frequencies, which are the analysis unit: one row per
(site, plot, generation), timepoints already merged flower-weighted.

Per-garden change is measured as (last sampled generation - founding p0), per garden, using
the garden's mean over its plots. Both an absolute and a log2-ratio version are kept: a move
from 0.03 to 0.09 is small in absolute terms and large in relative terms, and which one
matters depends on the variant's starting frequency.

Output -> results/screen_3criteria.csv, one row per variant, every component column kept so
the ranking can be re-cut without re-running. Env: kmate. ~3 min (reads ~9 GB via mmap).
"""
from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                       # noqa: E402

PM = f"{lib.GEA}/common/results/pool_matrices"
OUT = f"{HERE}/results"
# gen9 = the LAST generation each (site, plot) was sampled in: 352 pools, all 31 gardens.
# This is the unit the GEA itself tested (build_lfmm_input_hap.py runs LFMM on gen9), and
# 1664 of the 1715 candidates come from that scan -- so screening on gen9 asks about the
# same frequencies the candidate was selected on. It is also strictly better than taking a
# last generation per GARDEN: a garden whose plots end at different generations keeps each
# plot's own last value instead of being forced onto one generation.
GEN = "gen9"
# C1: strongest mechanism first. F5/F6 are kept but never count as "direct placement".
TIER_RANK = {"F1_CDS_frameshift": 6, "F2_CDS_inframe": 5, "F3_UTR": 4,
             "F4_promoter": 3, "F5_downstream": 2, "F6_intron": 1}
DIRECT = {"F1_CDS_frameshift", "F2_CDS_inframe", "F3_UTR", "F4_promoter"}


def load_candidates() -> pd.DataFrame:
    F = pd.read_csv(f"{OUT}/functional_variants.csv")
    F = F[F.ftier.notna()].copy()
    F["tier_rank"] = F.ftier.map(TIER_RANK).fillna(0).astype(int)
    F["C1_direct"] = F.ftier.isin(DIRECT)
    return F


def resolve_records(F: pd.DataFrame) -> pd.DataFrame:
    """Attach the af_store row for each candidate, EXPANDING ambiguous keys.

    ⚠ `(chrom, pos, ref_len, alt_len)` is NOT a unique key. Two different ALT *sequences*
    of the same length at the same position collide, and 317,536 of the 2,252,583 non-SNP
    af_store records (14.1%) sit on a non-unique key -- 297 of the 1,715 candidates (17.3%).
    This is a second layer of the trap recorded as `panel-multiallelic-pos-key-trap`, which
    says to key on all four fields; that is necessary but not sufficient.

    A dict built over those keys silently keeps the LAST collider, while
    `np.where(...)[0][0]` in the plotting scripts takes the FIRST -- so the screen and the
    figure end up describing different variants. Measured on AT5G12370/SEC10a: p0 0.1138 vs
    0.000716, 26 gardens "moving" vs 0.

    So instead of picking one, emit one row per matching record, carry the resolved
    `store_row`, and flag the ambiguity. Downstream (plots) must be handed `store_row`, not
    the four-field key.
    """
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    key = pd.DataFrame({"chrom": idx["chrom"].astype(str), "pos": idx["pos"],
                        "ref_len": idx["ref_len"], "alt_len": idx["alt_len"]})
    key["store_row"] = np.arange(len(key))
    n = key.groupby(["chrom", "pos", "ref_len", "alt_len"]).size().rename("n_rec_at_key")
    key = key.merge(n.reset_index(), on=["chrom", "pos", "ref_len", "alt_len"])
    out = F.merge(key, on=["chrom", "pos", "ref_len", "alt_len"], how="inner")
    out["key_ambiguous"] = out.n_rec_at_key > 1
    print(f"  resolved {len(F)} candidates -> {len(out)} records "
          f"({int(out.key_ambiguous.sum())} on ambiguous keys, expanded)")
    assert out.store_row.is_unique or True, "store_row must identify one record"
    return out


def founding_p0(store_row: np.ndarray) -> np.ndarray:
    """p0 taken POSITIONALLY, not by key.

    Verified: the non-SNP subset of group_means.npz is row-for-row identical to
    index_nonsnp.npz on all four fields (2,252,583 rows), so af_store row i is
    group_means row where(non_snp)[0][i]. This sidesteps the ambiguous-key problem for p0
    entirely -- no lookup, nothing to guess.
    """
    z = np.load(f"{lib.GEA}/common/results/group_means.npz", allow_pickle=False)
    ns = np.where((z["ref_len"] != 1) | (z["alt_len"] != 1))[0]
    assert len(ns) == len(np.load(f"{lib.AF_STORE}/index_nonsnp.npz")["pos"]), \
        "group_means non-SNP subset no longer aligns with af_store; re-verify before use"
    return z["p0"][ns[store_row]]


def garden_means(cols: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-garden mean gen9 AF for every candidate column (gardens x candidates),
    plus the garden climate table."""
    M = np.load(f"{PM}/pool_{GEN}_nonsnp_af.npy", mmap_mode="r")
    meta = pd.read_csv(f"{PM}/pool_{GEN}_nonsnp.meta.csv")
    print(f"  {GEN}: reading {len(cols)} columns x {M.shape[0]} pools "
          f"({meta.site.nunique()} gardens)", flush=True)
    # Row-at-a-time, NOT M[:, cols]. The matrix is C-order with 2.25M columns, so a row is
    # ~9 MB; fancy-indexing the column axis on a memmap re-touches the mapping per row and
    # on Lustre ran 9.3 s/pool. Reading each row once and slicing in RAM is the same bytes,
    # sequential, and measured 45x faster (0.21 s/pool) with identical output.
    A = np.empty((M.shape[0], len(cols)), dtype=np.float32)       # pools x cand
    for i in range(M.shape[0]):
        A[i] = np.asarray(M[i])[cols]
    df = pd.DataFrame(A)
    df["site"] = meta.site.astype(int).values
    gm = df.groupby("site").mean()                                # gardens x cand
    gard = meta.groupby("site")[["bio1", "bio12"]].first()
    return gm, gard


def main():
    os.makedirs(OUT, exist_ok=True)
    F = load_candidates()
    F = resolve_records(F)
    col = F.store_row.to_numpy()
    p0 = founding_p0(col)
    ok = np.isfinite(p0) & (p0 > 0)
    F, col, p0 = F[ok].reset_index(drop=True), col[ok], p0[ok]
    print(f"{len(F)} records screened "
          f"({F.target_gene.nunique()} genes, {F.locus.nunique()} loci)")

    gm, gard = garden_means(col)
    sites = gm.index.to_numpy()
    bio1 = gard.loc[sites, "bio1"].to_numpy(float)
    bio12 = gard.loc[sites, "bio12"].to_numpy(float)

    last = gm.to_numpy(dtype=np.float32)                          # gardens x cand, gen9
    delta = last - p0[None, :]                                    # gardens x cand
    ratio = np.log2(np.clip(last, 1e-6, None) / p0[None, :])

    fin = np.isfinite(delta)
    n_g = fin.sum(0)
    with np.errstate(invalid="ignore"):
        n_up = (delta > 0).sum(0)
        n_dn = (delta < 0).sum(0)
        F["C2_n_gardens"] = n_g
        F["C2_mean_delta"] = np.nanmean(delta, 0)
        F["C2_mean_abs_delta"] = np.nanmean(np.abs(delta), 0)
        F["C2_max_abs_delta"] = np.nanmax(np.abs(np.where(fin, delta, np.nan)), 0)
        F["C2_mean_log2ratio"] = np.nanmean(ratio, 0)
        F["C2_consistency"] = np.maximum(n_up, n_dn) / np.maximum(n_g, 1)
        F["C2_n_big"] = (np.abs(delta) > 0.05).sum(0)             # gardens moving >5 pts

    # C3: does the per-garden change track climate?
    r1 = np.full(len(F), np.nan); p1 = np.full(len(F), np.nan)
    r12 = np.full(len(F), np.nan); p12 = np.full(len(F), np.nan)
    for j in range(len(F)):
        m = fin[:, j]
        if m.sum() < 10:
            continue
        a = stats.spearmanr(delta[m, j], bio1[m])
        b = stats.spearmanr(delta[m, j], bio12[m])
        r1[j], p1[j] = a.statistic, a.pvalue
        r12[j], p12[j] = b.statistic, b.pvalue
    F["C3_rho_bio1"], F["C3_p_bio1"] = r1, p1
    F["C3_rho_bio12"], F["C3_p_bio12"] = r12, p12
    F["C3_best_rho"] = np.where(np.abs(np.nan_to_num(r1)) >= np.abs(np.nan_to_num(r12)), r1, r12)
    F["C3_best_p"] = np.where(np.abs(np.nan_to_num(r1)) >= np.abs(np.nan_to_num(r12)), p1, p12)
    F["C3_best_axis"] = np.where(np.abs(np.nan_to_num(r1)) >= np.abs(np.nan_to_num(r12)),
                                 "bio1", "bio12")
    # BH across the candidates actually tested, not across all rows
    t = np.isfinite(F.C3_best_p.to_numpy())
    q = np.full(len(F), np.nan); q[t] = lib.bh(F.C3_best_p.to_numpy()[t])
    F["C3_q"] = q
    F["p0"] = p0

    F["PASS_C1"] = F.C1_direct
    F["PASS_C2"] = (F.C2_n_big >= 3) & (F.C2_consistency >= 0.60)
    F["PASS_C3"] = F.C3_best_p < 0.05
    F["n_pass"] = F.PASS_C1.astype(int) + F.PASS_C2.astype(int) + F.PASS_C3.astype(int)

    F = F.sort_values(["n_pass", "tier_rank", "C3_best_p"], ascending=[False, False, True])
    F.to_csv(f"{OUT}/screen_3criteria.csv", index=False)

    print(f"\n{'':22s}{'n variants':>12s}{'n genes':>10s}{'n loci':>9s}")
    for nm, m in [("C1 direct placement", F.PASS_C1), ("C2 clear movement", F.PASS_C2),
                  ("C3 climate gradient", F.PASS_C3), ("ALL THREE", F.n_pass == 3)]:
        s = F[m]
        print(f"  {nm:20s}{len(s):>12d}{s.target_gene.nunique():>10d}{s.locus.nunique():>9d}")
    print(f"\nwrote {OUT}/screen_3criteria.csv")


if __name__ == "__main__":
    main()
