#!/usr/bin/env python
"""Per-founder SV load under different record filters, to reconcile `sv_origin_partition` with
`sv_founder_confound`.

The founder-level result ("cold-origin founders carry more inserted kb", rho = -0.54) is computed
on COMMON, well-called insertions only (MAC >= 12 and call rate >= 0.9, `_founder_sv_content.py`).
The origin partition uses every SV record. This script computes each founder's carried kb and
count under several filters, and split by carrier-count stratum, so the two can be compared.

Record-level filters (match `_founder_sv_content.py`): common_wellcalled, common, wellcalled, rare,
rare_wellcalled, all. LOCUS-level filters and strata (`*_min2`, `*_strata`): records of one SV type at one position are one
locus (alternative alleles of the same event); `*_min2` and the strata use the locus carrier count
(union of carriers over its records). A founder carries at most one record per position, so
per-founder totals are unaffected by the grouping. RECORD-level versions (`*_rec2`,
`*_strata_rec`) use the allele's own carrier count; they drop private alleles of shared loci,
which are almost all cactus-assembly sequence variants. Founder covariates (SNP divergence, SV call
rate) are copied from founder_sv_content.npz for the partial correlations in the notebook.

Env: kmate. Writes results/sv_adaptive/sv_origin_load_by_filter.npz.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
SV_BP, MIN_MAC, CALL_MIN = 50, 12, 0.9
FILTERS = ["common_wellcalled", "common", "wellcalled", "rare", "rare_wellcalled", "all",
           "rare_min2", "rare_min2_wellcalled", "all_min2",   # *_min2: LOCUS has >= 2 carriers (rare: 2-11)
           "rare_rec2", "rare_rec2_wellcalled", "all_rec2"]   # *_rec2: RECORD (allele) has >= 2 carriers (rare: 2-11)
STRATA = [(1, 1), (2, 3), (4, 7), (8, 20), (21, 77), (78, 10 ** 6)]


def fmask(name, a, c, n, nl):
    # a, c: record carriers / called; nl: carriers of the record's locus
    common = (a >= MIN_MAC) & (a <= n - MIN_MAC)
    rare = (a >= 1) & (a < MIN_MAC)
    rare2 = (nl >= 2) & (nl < MIN_MAC)
    well = c >= CALL_MIN * n
    return {"common_wellcalled": common & well, "common": common, "wellcalled": (a >= 1) & well,
            "rare": rare, "rare_wellcalled": rare & well, "all": a >= 1,
            "rare_min2": rare2, "rare_min2_wellcalled": rare2 & well, "all_min2": nl >= 2,
            "rare_rec2": (a >= 2) & (a < MIN_MAC), "rare_rec2_wellcalled": (a >= 2) & (a < MIN_MAC) & well,
            "all_rec2": a >= 2}[name]


def main():
    F = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    n = len(F)
    acc = {}
    for t in ("ins", "del"):
        for f in FILTERS:
            acc[f"kb_{t}_{f}"] = np.zeros(n); acc[f"n_{t}_{f}"] = np.zeros(n); acc[f"nrec_{t}_{f}"] = 0
        acc[f"kb_{t}_strata"] = np.zeros((len(STRATA), n)); acc[f"n_{t}_strata"] = np.zeros((len(STRATA), n)); acc[f"kb_{t}_strata_rec"] = np.zeros((len(STRATA), n)); acc[f"n_{t}_strata_rec"] = np.zeros((len(STRATA), n)); acc[f"nrec_{t}_strata"] = np.zeros(len(STRATA), int)
    for ch in lib.CHROMS:
        cl = ch.lower(); base = f"{lib.PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        assert np.array_equal(meta["founders"].astype("U6"), F.founder.to_numpy()), f"{ch}: founder order differs"
        rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsr(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        a = np.asarray(vp.sum(0)).ravel(); c = np.asarray(vc.sum(0)).ravel()
        pos = meta["pos"].astype(np.int64)
        for t, sel in (("ins", (al - rl) > SV_BP), ("del", (rl - al) > SV_BP)):
            size = np.abs(al - rl).astype(float)
            # locus carrier count per record (records of this type at the same position)
            cols = np.where(sel & (a >= 1))[0]
            upos, li = np.unique(pos[cols], return_inverse=True); li = li.ravel()
            Mloc = sp.csr_matrix((np.ones(cols.size), (li, np.arange(cols.size))), shape=(upos.size, cols.size))
            nloc = np.asarray(((Mloc @ vp[:, cols].T.tocsc()) > 0).sum(1)).ravel()
            nl = np.zeros(a.size); nl[cols] = nloc[li]
            for f in FILTERS:
                m = sel & fmask(f, a, c, n, nl)
                acc[f"kb_{t}_{f}"] += vp @ np.where(m, size, 0) / 1e3
                acc[f"n_{t}_{f}"] += vp @ m.astype(float)
                acc[f"nrec_{t}_{f}"] += int(m.sum())
            for i, (lo, hi) in enumerate(STRATA):
                m = sel & (nl >= lo) & (nl <= hi)
                acc[f"kb_{t}_strata"][i] += vp @ np.where(m, size, 0) / 1e3
                acc[f"n_{t}_strata"][i] += vp @ m.astype(float)
                m = sel & (a >= lo) & (a <= hi)                       # same, by RECORD carrier count
                acc[f"kb_{t}_strata_rec"][i] += vp @ np.where(m, size, 0) / 1e3
                acc[f"n_{t}_strata_rec"][i] += vp @ m.astype(float)
                acc[f"nrec_{t}_strata"][i] += int(m.sum())
        print(f"[{ch}] done", flush=True)
    S = np.load(f"{OUT}/founder_sv_content.npz", allow_pickle=True)
    assert np.array_equal(S["founders"].astype(str), F.founder.to_numpy().astype(str))
    for k in ("n_snp_all", "n_snp", "call_rate_sv", "call_rate_all"):
        acc[f"cov_{k}"] = S[k].astype(float)
    np.savez_compressed(f"{OUT}/sv_origin_load_by_filter.npz", founders=F.founder.to_numpy().astype(str), bio1=F.bio1.to_numpy().astype(float),
                        group=F.group.to_numpy().astype(str), filters=np.array(FILTERS, dtype=str), strata=np.array(STRATA), **acc)
    print(f"[wrote] {OUT}/sv_origin_load_by_filter.npz")


if __name__ == "__main__":
    main()
