#!/usr/bin/env python
"""GWAS re-screen, v2: represent every GWAS-involved gene by the allele that was significant.

`screen_gwas_rescreen.csv` (build_ownaxis_shortlist.py) picked each gene's record by the
screen criteria, and the GWAS garden counts it carries were joined by POSITION through
`build_gwas_pool.py`, which assigns the largest panel record at a multiallelic position.
`resolve_gwas_alleles.py` showed that allele is wrong at 21 of the pool's 33 multiallelic
markers, so 10 of the 45 genes were drawn with a non-significant allele (CYP71B4, CML25,
EPFL5 showed 0/30 gardens), and 3 GEA+GWAS genes were drawn with their GEA record.

Here, per gene:
  * GWAS positions = positions of the gene's screen records that carry a significant allele
    in `gwas_hits_allele_resolved.csv`
  * record = the af_store row with that exact (chrom, pos, ref_len, alt_len) AND the same
    rank among colliding same-length alleles as the significant panel record (`key_rank`)
  * representative = the allele significant in the most gardens; ties -> higher GWAS nlp
  * movement: the same relative criterion as v1 (C2rel), on gen9 garden means
  * climate: Spearman of (gen9 garden mean - p0) on the axis v1 used, as context only
  * GWAS direction: sign of the GEMMA Z in the significant gardens (+ = carriers have
    higher s) next to the direction the pool AF moved in those same gardens

Output -> results/screen_gwas_rescreen_v2.csv (one row per gene; `rep_changed` marks the
genes whose drawn allele differs from v1).

env: kmate. Compute node. ~3-5 min (gen9 pool matrix read).
"""
from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
import lib                                                       # noqa: E402
import axis_clusters as ac                                       # noqa: E402
from screen_3criteria import founding_p0                        # noqa: E402
from build_ownaxis_shortlist import garden_means                # noqa: E402

OUT = f"{HERE}/results"


def tier(old: str, rl: int, al: int) -> str:
    """Placement is positional; only the CDS frame call depends on the allele size."""
    if isinstance(old, str) and old.startswith(("F1_", "F2_")):
        return "F1_CDS_frameshift" if abs(al - rl) % 3 else "F2_CDS_inframe"
    return old


def main():
    V1 = pd.read_csv(f"{OUT}/screen_gwas_rescreen.csv")
    H = pd.read_csv(f"{OUT}/gwas_hits_allele_resolved.csv")
    H["chrom"] = H.chrom.str.replace("chr", "Chr")
    # an allele can be significant in both the nonsnp and the sv scan: keep the wider one
    H = (H.sort_values("n_gardens", ascending=False)
          .drop_duplicates(["chrom", "pos", "ref_len", "alt_len"]))
    S = pd.read_csv(f"{OUT}/screen_3criteria.csv",
                    usecols=["target_gene", "chrom", "pos", "ftier", "region"])
    S = S[S.target_gene.isin(V1.target_gene)].drop_duplicates(["target_gene", "chrom", "pos"])
    C = S.merge(H, on=["chrom", "pos"])                      # gene x significant allele

    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    key = pd.DataFrame({"chrom": idx["chrom"].astype(str), "pos": idx["pos"],
                        "ref_len": idx["ref_len"], "alt_len": idx["alt_len"]})
    key["store_row"] = np.arange(len(key))
    key = key[key.pos.isin(C.pos)]
    key["key_rank"] = key.groupby(["chrom", "pos", "ref_len", "alt_len"]).cumcount()
    n_at = key.groupby(["chrom", "pos", "ref_len", "alt_len"]).store_row.transform("count")
    key["key_ambiguous"] = n_at > 1
    # same-length ALT alleles collide on the 4-field key: take the store row with the same
    # rank as the significant panel record (store and panel orders agree), never the one
    # that moved most -- that picked CRK18's non-significant allele in the first version
    C = C.merge(key, on=["chrom", "pos", "ref_len", "alt_len", "key_rank"], how="left")
    miss = C[C.store_row.isna()]
    if len(miss):
        print(f"no af_store record for {len(miss)} significant alleles:\n"
              f"{miss[['target_gene', 'chrom', 'pos', 'ref_len', 'alt_len']]}")
    C = C.dropna(subset=["store_row"]).copy()
    C["store_row"] = C.store_row.astype(int)
    C["p0"] = founding_p0(C.store_row.to_numpy())

    gm = garden_means(np.unique(C.store_row.to_numpy()))
    sites = gm.index.to_numpy()
    A = ac.site_climate()
    A.index = sorted(pd.read_csv(ac.POOLMETA).site.astype(int).unique())
    Ax = A.loc[sites]
    axis_of = V1.set_index("target_gene").C3_axis_used

    rows = []
    for r in C.itertuples():
        last = gm[r.store_row].to_numpy(float)
        delta = last - r.p0
        ratio = np.log2(np.clip(last, 1e-6, None) / max(r.p0, 1e-6))
        up = (ratio >= 1) & (delta >= 0.02)
        dn = (ratio <= -1) & (delta <= -0.02)
        sig = [int(g) for g in str(r.gardens).split(",")]
        zs = [float(z) for z in str(r.Z_sig_gardens).split(",")]
        ins = np.isin(sites, sig)
        d_sig = delta[ins]
        ax = axis_of.get(r.target_gene, "bio1")
        m = np.isfinite(delta)
        c3 = stats.spearmanr(delta[m], Ax[ax].to_numpy(float)[m])
        rows.append({
            "target_gene": r.target_gene, "chrom": r.chrom, "pos": r.pos,
            "ref_len": r.ref_len, "alt_len": r.alt_len, "size": abs(r.alt_len - r.ref_len),
            "gwas_cls": r.cls, "ftier": tier(r.ftier, r.ref_len, r.alt_len),
            "region": r.region, "p0": r.p0, "store_row": r.store_row,
            "key_ambiguous": r.key_ambiguous,
            "gwas_n_gardens": r.n_gardens, "gwas_gardens": r.gardens,
            "gwas_best_nlp": r.best_nlp, "gwas_Z_sign": "+" if np.mean(zs) > 0 else "-",
            "sig_gardens_af_delta": ",".join(f"{v:+.3f}" for v in d_sig),
            "sig_gardens_n_up": int((up & ins).sum()), "sig_gardens_n_dn": int((dn & ins).sum()),
            "C2rel_n_up": int(up.sum()), "C2rel_n_dn": int(dn.sum()),
            "C2rel_max_last": float(np.nanmax(last)),
            "PASS_C2rel": bool(max(up.sum(), dn.sum()) >= 3
                               and max(up.sum(), dn.sum()) / max(up.sum() + dn.sum(), 1) >= 0.6),
            "C3_axis_used": ax, "C3_rho_used": c3.statistic, "C3_p_used": c3.pvalue,
            "PASS_gwas2": r.n_gardens >= 2})
    R = pd.DataFrame(rows)
    # ties broken by GWAS strength, never by how much the allele moved (that would select
    # the representative on the outcome being judged)
    R = (R.sort_values(["gwas_n_gardens", "gwas_best_nlp"], ascending=False)
          .drop_duplicates("target_gene"))

    keep = V1.set_index("target_gene")
    R["symbol"] = R.target_gene.map(keep.symbol)
    R["found_by"] = R.target_gene.map(keep.found_by)
    R["protein"] = R.target_gene.map(keep.protein)
    R["r_garden4_set_v1"] = R.target_gene.map(keep.r_garden4_set)
    v1key = keep.chrom + ":" + keep.pos.astype(str) + ":" + keep.ref_len.astype(str) + ":" \
        + keep.alt_len.astype(str)
    R["rep_changed"] = (R.chrom + ":" + R.pos.astype(str) + ":" + R.ref_len.astype(str) + ":"
                        + R.alt_len.astype(str)) != R.target_gene.map(v1key)
    lost = sorted(set(V1.target_gene) - set(R.target_gene))
    first = ["target_gene", "symbol", "found_by", "rep_changed"]
    R = R[first + [c for c in R.columns if c not in first]]
    R = R.sort_values(["PASS_gwas2", "gwas_n_gardens"], ascending=False)
    R.to_csv(f"{OUT}/screen_gwas_rescreen_v2.csv", index=False)
    print(f"wrote screen_gwas_rescreen_v2.csv: {len(R)} genes, rep changed for "
          f"{int(R.rep_changed.sum())}, PASS_gwas2 {int(R.PASS_gwas2.sum())} "
          f"(v1: {int(V1.PASS_gwas2.sum())}); genes with no resolvable allele: {lost}")
    pd.set_option("display.width", 250)
    print(R[["symbol", "rep_changed", "ref_len", "alt_len", "ftier", "p0", "gwas_n_gardens",
             "gwas_gardens", "gwas_Z_sign", "sig_gardens_af_delta", "C2rel_n_up",
             "C2rel_n_dn", "C3_axis_used", "C3_rho_used", "key_ambiguous"]]
          .round(3).to_string(index=False))


if __name__ == "__main__":
    main()
