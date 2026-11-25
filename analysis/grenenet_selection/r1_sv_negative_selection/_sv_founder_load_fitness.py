#!/usr/bin/env python
"""Founder-level test: do founders carrying more cold-specific insertions lose frequency in GrENE-Net
beyond what their origin climate and sequencing route predict?

Per-variant selection is not identifiable here: pool frequencies are a global-mode projection onto
founder genotypes, so a variant's trajectory is fixed by its carrier set. What IS identifiable is each
founder's own frequency change. So the purging question is asked of founders:
  fitness_f ~ load_f + origin bio1_f + cactus_f   (+ restricted-indel count as the matched control)

FITNESS (from r3_persite_gwas/ecotype_fitness.py, founder_site_dh.npz: DH_rel [31 sites x 231], H0):
  logfit  = mean over gardens of log((h0 + dh + eps) / (h0 + eps)), eps = 1e-5 -- scale-free; founders
            start at frequencies spanning ~75-fold, so absolute dh is dominated by common founders
  fdecl   = fraction of gardens in which the founder's frequency falls (dh < 0)
LOAD, from the record-level (fixed) allele set `sv_purging_classify_variants_rec2.csv`:
  number of cold-specific / warm-specific / shared insertion alleles carried, and the same for indels.

Env: kmate. Writes results/sv_adaptive/sv_founder_load_fitness.csv.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
EF = f"{lib.GEA}/r3_persite_gwas/results/ecotype_fitness/founder_site_dh.npz"
EPS = 1e-5


def main():
    F = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    Z = np.load(EF, allow_pickle=True)
    assert np.array_equal(Z["founders"].astype(str), F.founder.to_numpy()), "founder order differs"
    DH = Z["DH_rel"]; H0 = Z["H0"]                       # [sites x founders], [founders]
    F["logfit"] = np.nanmean(np.log((H0[None, :] + DH + EPS) / (H0[None, :] + EPS)), axis=0)
    F["fdecl"] = np.nanmean(DH < 0, axis=0)
    F["h0"] = H0

    V = pd.read_csv(f"{OUT}/sv_purging_classify_variants_rec2.csv")
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz"); ch = idx["chrom"].astype("U5")
    off, o = {}, 0
    for C in lib.CHROMS:
        off[C] = o; o += int((ch == C).sum())
    groups = [(c, g) for c in ("insertion", "indel", "deletion") for g in ("cold-specific", "warm-specific", "shared")]
    for c, g in groups:
        F[f"n_{c[:3]}_{g.split('-')[0]}"] = 0
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        prl = meta["ref_len"].astype(np.int64); pal = meta["alt_len"].astype(np.int64)
        nonsnp = np.where(~((prl == 1) & (pal == 1)))[0]
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
        for c, g in groups:
            d = V[(V.cls == c) & (V.origin == g) & (V.chrom == C)]
            if len(d) == 0:
                continue
            cols = nonsnp[d.afrow.to_numpy() - off[C]]
            F[f"n_{c[:3]}_{g.split('-')[0]}"] += np.asarray((vp[:, cols] > 0).sum(1)).ravel()
        print(f"[{C}] done", flush=True)
    F.to_csv(f"{OUT}/sv_founder_load_fitness.csv", index=False)
    print(f"[wrote] {OUT}/sv_founder_load_fitness.csv")


if __name__ == "__main__":
    main()
