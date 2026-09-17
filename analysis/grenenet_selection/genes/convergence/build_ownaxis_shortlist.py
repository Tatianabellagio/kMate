#!/usr/bin/env python
"""Build the candidate lists that replace the bio1/bio12-only cut (user decision 2026-09-16).

Two outputs, both one row per gene, both carrying `store_row` so the plots draw the same
record the table describes:

  results/screen_top_loci_ownaxis.csv
      Every gene passing C1 (direct placement) + C2 (movement) + C3 judged on the axis the
      GEA actually found it on (`screen_c3_all_axes_full.csv`), that is NOT already in the
      first visual review (`screen_visual_review.csv`). Representative record = the passing
      record with the smallest own-axis p.

  results/screen_gwas_rescreen.csv
      Every GWAS-involved gene (found_by GWAS or GEA+GWAS) re-screened with a movement
      criterion suited to rare alleles. GWAS hits are MAC >= 5 by construction, so
      `C2_n_big` (>= 3 gardens moving > 5 points) cannot be met by an allele at p0 0.01
      even under strong selection -- CRK18 fails it for exactly that reason. Here:
        C2rel  >= 3 gardens with log2(last/p0) >= 1 AND (last - p0) >= 0.02, consistency >= 0.6
        C3     GEA+GWAS genes: own GEA axis.  GWAS-only genes have no a-priori axis, so the
               best of the 22 is reported and flagged `C3_exploratory` -- a 22-way search,
               not a confirmation.

Both tables get `r_garden4_set`: Spearman correlation of the gene's per-garden gen9 mean
profile with the mean profile of the 20 records flagged as the garden-4 founder sweep in the
first review. Founder-sweep alleles co-vary at r ~0.9 across chromosomes; a new candidate at
r > 0.6 is riding that haplotype, not a locus signal.

Env: kmate. Compute node. ~5-10 min (one Lustre pass over the gen9 pool matrix).
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
import axis_clusters as ac                                       # noqa: E402

PM = f"{lib.GEA}/common/results/pool_matrices"
OUT = f"{HERE}/results"
GEN = "gen9"


def garden_means(cols: np.ndarray) -> pd.DataFrame:
    M = np.load(f"{PM}/pool_{GEN}_nonsnp_af.npy", mmap_mode="r")
    meta = pd.read_csv(f"{PM}/pool_{GEN}_nonsnp.meta.csv")
    A = np.empty((M.shape[0], len(cols)), dtype=np.float32)
    for i in range(M.shape[0]):
        A[i] = np.asarray(M[i])[cols]
    df = pd.DataFrame(A, columns=cols)
    df["site"] = meta.site.astype(int).values
    return df.groupby("site").mean()


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--gea-only", action="store_true",
                    help="write only screen_top_loci_ownaxis.csv; leave the (hand-annotated) "
                         "screen_gwas_rescreen.csv alone -- superseded by fix_gwas_reps.py")
    a = ap.parse_args()
    R = pd.read_csv(f"{OUT}/screen_c3_all_axes_full.csv")
    S = pd.read_csv(f"{OUT}/screen_3criteria.csv")
    M = pd.read_csv(f"{OUT}/master_candidate_genes.csv").set_index("gene")
    V = pd.read_csv(f"{OUT}/screen_visual_review.csv", dtype=str).fillna("")
    reviewed = set(V.target_gene)

    # per-record extras from the 3-criteria table, keyed on store_row (unique per record)
    extra = S.set_index("store_row")[["ref_len", "alt_len", "region", "size", "vclass",
                                       "in_gwas", "gwas_n_gardens", "gwas_gardens",
                                       "gwas_mac", "locus", "locus_n_genes",
                                       "snp_cosig_2kb", "best_r2_snp", "te_overlap"]]
    R = R.join(extra, on="store_row")
    R["found_by"] = R.target_gene.map(M.found_by)
    # same-length colliders on the 4-field key: keep only records the GEA actually flagged
    # (audit_gea_colliders.py, 2026-09-17). 805 of 1,111 expanded records were not.
    G_ = pd.read_csv(f"{OUT}/screen_record_gea_match.csv").drop_duplicates("store_row")
    R["gea_sig_record"] = R.store_row.map(G_.set_index("store_row").gea_sig_record)
    R["symbol"] = R.target_gene.map(M.symbol)

    # ---------- 1. own-axis shortlist: new genes only ----------
    ok = R.PASS_C1 & R.PASS_C2 & R.PASS_C3own & (R.gea_sig_record != False)
    P = R[ok & ~R.target_gene.isin(reviewed)].copy()
    P["abs_rho"] = P.C3own_rho.abs()
    P = P.sort_values(["C3own_p", "abs_rho"], ascending=[True, False])
    new = P.drop_duplicates("target_gene").copy()
    print(f"own-axis shortlist: {len(new)} new genes "
          f"({new.found_by.value_counts().to_dict()}); axes {new.gea_best_axis.value_counts().to_dict()}")

    if a.gea_only:
        G = R.iloc[:0]
    # ---------- 2. GWAS re-screen ----------
    G = R[R.found_by.isin(["GWAS", "GEA+GWAS"])].copy()
    print(f"GWAS-involved records: {len(G)} ({G.target_gene.nunique()} genes)")

    # ---------- garden means for everything we need ----------
    g4_rows = V[V.notes.str.contains("garden-4", na=False)].store_row.astype(int).to_numpy()
    cols = np.unique(np.concatenate([new.store_row.to_numpy(), G.store_row.to_numpy(), g4_rows]))
    print(f"reading {len(cols)} columns of the gen9 pool matrix ...", flush=True)
    gm = garden_means(cols)
    sites = gm.index.to_numpy()
    g4_profile = gm[g4_rows].mean(axis=1).to_numpy()

    A = ac.site_climate()
    A.index = sorted(pd.read_csv(ac.POOLMETA).site.astype(int).unique())
    Ax = A.loc[sites]

    def r_g4(row):
        return float(stats.spearmanr(gm[row].to_numpy(), g4_profile).statistic)

    new["r_garden4_set"] = [round(r_g4(r), 2) for r in new.store_row]
    new["score"] = new.abs_rho * new.C2_consistency
    new = new.sort_values("score", ascending=False)
    keep = ["target_gene", "symbol", "found_by", "chrom", "pos", "ref_len", "alt_len", "size",
            "vclass", "ftier", "region", "p0", "C2_n_big", "C2_consistency", "gea_best_axis",
            "gea_sig_axes", "C3own_axis", "C3own_rho", "C3own_p", "C3all_best_axis",
            "C3all_best_rho", "C3all_same_cluster", "rho_bio1", "rho_bio12", "r_garden4_set",
            "score", "store_row", "key_ambiguous", "locus", "locus_n_genes", "snp_cosig_2kb",
            "best_r2_snp", "te_overlap", "in_gwas", "gwas_n_gardens"]
    new[keep].to_csv(f"{OUT}/screen_top_loci_ownaxis.csv", index=False)
    print(f"wrote {OUT}/screen_top_loci_ownaxis.csv ({len(new)} rows); "
          f"r_garden4>0.6: {(new.r_garden4_set > 0.6).sum()}")

    if a.gea_only:
        return
    # GWAS: relative movement per record
    last = gm[G.store_row.to_numpy()].to_numpy().T                   # records x gardens
    p0 = G.p0.to_numpy(np.float32)[:, None]
    delta = last - p0
    ratio = np.log2(np.clip(last, 1e-6, None) / np.clip(p0, 1e-6, None))
    up = (ratio >= 1) & (delta >= 0.02)
    dn = (ratio <= -1) & (delta <= -0.02)
    n_up, n_dn = up.sum(1), dn.sum(1)
    G["C2rel_n_up"], G["C2rel_n_dn"] = n_up, n_dn
    G["C2rel_max_last"] = np.nanmax(last, axis=1)
    G["C2rel_consistency"] = np.maximum(n_up, n_dn) / np.maximum(n_up + n_dn, 1)
    G["PASS_C2rel"] = (np.maximum(n_up, n_dn) >= 3) & (G.C2rel_consistency >= 0.6)
    both = G.found_by == "GEA+GWAS"
    G["C3_axis_used"] = np.where(both, G.gea_best_axis, G.C3all_best_axis)
    G["C3_rho_used"] = np.where(both, G.C3own_rho, G.C3all_best_rho)
    G["C3_p_used"] = np.where(both, G.C3own_p, G.C3all_best_p)
    G["C3_exploratory"] = ~both
    G["PASS_C3used"] = G.C3_p_used < 0.05
    G["PASS_new"] = G.PASS_C1 & G.PASS_C2rel & G.PASS_C3used
    G["r_garden4_set"] = [round(r_g4(r), 2) for r in G.store_row]
    G["abs_rho"] = G.C3_rho_used.abs()
    G = G.sort_values(["PASS_new", "PASS_C1", "PASS_C2rel", "abs_rho"],
                      ascending=[False, False, False, False])
    rep = G.drop_duplicates("target_gene")
    keepg = ["target_gene", "symbol", "found_by", "chrom", "pos", "ref_len", "alt_len", "size",
             "vclass", "ftier", "region", "p0", "gwas_n_gardens", "gwas_gardens", "gwas_mac",
             "PASS_C1", "PASS_C2", "C2rel_n_up", "C2rel_n_dn", "C2rel_max_last",
             "C2rel_consistency", "PASS_C2rel", "C3_axis_used", "C3_rho_used", "C3_p_used",
             "C3_exploratory", "PASS_C3used", "PASS_new", "r_garden4_set", "store_row",
             "key_ambiguous", "locus", "snp_cosig_2kb", "best_r2_snp"]
    rep[keepg].to_csv(f"{OUT}/screen_gwas_rescreen.csv", index=False)
    print(f"wrote {OUT}/screen_gwas_rescreen.csv ({len(rep)} genes); "
          f"pass C1+C2rel+C3: {int(rep.PASS_new.sum())}, "
          f"C1+C2rel only: {int((rep.PASS_C1 & rep.PASS_C2rel).sum())}, "
          f"C2rel any: {int(rep.PASS_C2rel.sum())}")
    print(rep[rep.PASS_C2rel][["target_gene", "symbol", "found_by", "ftier", "p0",
                               "C2rel_n_up", "C2rel_n_dn", "C2rel_max_last", "C3_axis_used",
                               "C3_rho_used", "C3_exploratory", "PASS_new",
                               "r_garden4_set"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
