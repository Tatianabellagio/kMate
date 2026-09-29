#!/usr/bin/env python
"""Per-garden GWAS markers, judged on their own terms -- not on a climate gradient.

The evidence matrix's climate lines reward a gradient across gardens, which a per-garden
GWAS hit is not expected to have: it is selection on the founder coefficient `s` in
particular gardens, possibly with no climate pattern at all (round-2 brief: "rare alleles
cannot give a clean climate gradient, and that is fine"). The brief's own criterion:
significant gardens, and whether the allele actually moves in them, consistently.

Per allele-resolved marker (store row by key_rank, as resolve_gwas_alleles.py):
  n_sig_gardens     Bonferroni gardens (allele-resolved)
  delta_sig         mean change from p0 over the pools of the significant gardens
  consistency       fraction of those pools moving in the same direction as delta_sig
  rises_over_time   Spearman rho of change vs generation within the significant gardens
  specificity_p     Mann-Whitney, pools of significant gardens vs all other pools
  z_agrees          sign of GEMMA's best Z agrees with the direction of delta_sig
  grade             STRONG >=3 gardens + consistent (>=0.75) + specific (p<0.01);
                    GOOD 2 gardens consistent+specific, or >=3 with one of the two;
                    MODERATE 2 gardens with one of the two; WEAK otherwise
Joined to the evidence matrix's mechanism / chromatin / expression / flags.
Writes results/gwas_evaluation.csv.
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", ".."))); sys.path.insert(0, HERE)
import lib, axis_clusters as ac                                  # noqa: E402

OUT = f"{HERE}/results"
K4 = ["chrom", "pos", "ref_len", "alt_len"]


def main():
    G = pd.read_csv(f"{OUT}/gwas_hits_allele_resolved.csv")
    G["chrom"] = G.chrom.str.replace("^chr", "Chr", regex=True)
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    K = pd.DataFrame({"chrom": idx["chrom"].astype(str), "pos": idx["pos"],
                      "ref_len": idx["ref_len"], "alt_len": idx["alt_len"]})
    K["row"] = np.arange(len(K))
    Kg = K.groupby(K4).row.apply(lambda x: np.sort(x.values)).to_dict()
    G["store_row"] = [int(Kg[(r.chrom, r.pos, r.ref_len, r.alt_len)][r.key_rank])
                      if (r.chrom, r.pos, r.ref_len, r.alt_len) in Kg else -1
                      for r in G.itertuples()]
    M = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy", mmap_mode="r")
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    PM = pd.read_csv(ac.POOLMETA)
    site, gen = PM.site.values.astype(int), PM.generation.values
    rows = []
    for g in G.itertuples():
        if g.store_row < 0:
            continue
        d = np.asarray(M[:, g.store_row], float) - p0[g.store_row]
        sig = {int(x) for x in str(g.gardens).split(",") if x.strip()}
        ins = np.isin(site, list(sig)); ok = np.isfinite(d)
        ds, do = d[ins & ok], d[~ins & ok]
        delta = float(np.mean(ds)) if len(ds) else np.nan
        cons = float(np.mean(np.sign(ds) == np.sign(delta))) if len(ds) else np.nan
        rho = spearmanr(gen[ins & ok], ds)[0] if len(ds) > 3 else np.nan
        spec = mannwhitneyu(ds, do).pvalue if len(ds) > 2 else np.nan
        rows.append(dict(store_row=g.store_row, chrom=g.chrom, pos=g.pos, ref_len=g.ref_len,
                         alt_len=g.alt_len, cls=g.cls, mac=g.mac, n_sig_gardens=g.n_gardens,
                         gardens=g.gardens, best_nlp=g.best_nlp, best_Z=g.best_Z, p0=p0[g.store_row],
                         delta_sig=delta, delta_other=float(np.mean(do)), consistency=cons,
                         rises_over_time=rho, specificity_p=spec,
                         z_agrees=bool(np.sign(g.best_Z) == np.sign(delta))))
    E = pd.DataFrame(rows)
    # DIRECTIONAL, not generic: rare alleles (MAC 5-9, p0 1-3%) drift down in every garden,
    # so "consistent decline" is the null expectation. What the GWAS claims is a direction
    # (sign of GEMMA's Z on s); the allele must move THAT way in its significant gardens,
    # consistently across their pools, and more than in the other gardens.
    E["direction"] = np.sign(E.best_Z)
    E["excess"] = (E.delta_sig - E.delta_other) * E.direction
    cons_dir = []
    for g, r in zip(G[G.store_row >= 0].itertuples(), E.itertuples()):
        d = np.asarray(M[:, r.store_row], float) - p0[r.store_row]
        sig = {int(x) for x in str(r.gardens).split(",") if x.strip()}
        ds = d[np.isin(site, list(sig)) & np.isfinite(d)]
        cons_dir.append(float(np.mean(np.sign(ds) == r.direction)) if len(ds) else np.nan)
    E["consistency_dir"] = cons_dir
    # consistency is REQUIRED for any grade above WEAK: a Mann-Whitney over hundreds of
    # pools is significant when a handful of pools spike, which is the brief's "single-pool
    # spike" pattern, not a garden-wide response
    good = (E.consistency_dir >= 0.6) & (E.excess > 0)
    spec = (E.specificity_p < 0.01) & (E.excess > 0) & (E.consistency_dir >= 0.5)
    E["grade"] = np.select(
        [(E.n_sig_gardens >= 3) & good & spec,
         ((E.n_sig_gardens == 2) & good & spec) | ((E.n_sig_gardens >= 3) & good),
         (E.n_sig_gardens == 2) & good],
        ["STRONG", "GOOD", "MODERATE"], "WEAK")
    X = pd.read_csv(f"{OUT}/evidence_matrix.csv")
    keep = ["store_row", "target_gene", "sym", "mode", "dist_tss", "mech_grade", "atac_bp",
            "expr_fold", "q_expr", "p_expr_emp", "L_expression", "local_rank", "snp_rivals",
            "carriers", "called", "af_shared_vector", "call_rate_lt90", "gene_silent"]
    E = E.merge(X[[c for c in keep if c in X.columns]], on="store_row", how="left")
    E.to_csv(f"{OUT}/gwas_evaluation.csv", index=False)
    print(E.grade.value_counts().to_string())


if __name__ == "__main__":
    main()
