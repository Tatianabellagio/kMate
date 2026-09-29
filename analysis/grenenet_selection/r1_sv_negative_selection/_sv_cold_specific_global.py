#!/usr/bin/env python
"""Do COLD-SPECIFIC insertions decline in EVERY garden, or only the warm ones?

The distinction the temporal result turns on:
  * true purifying selection -> a deleterious insertion loses ground everywhere, cold gardens too;
  * local adaptation / hitchhiking -> it falls where its carrier lineages lose and rises where they win.

`_sv_temporal_climate_matched.py` already wrote, per garden, each insertion's residual against SNPs
matched on p0 AND carrier-origin climate. So the comparison SNPs for a cold-carried insertion are
themselves cold-carried, and any residual decline is specific to the insertion rather than to its
lineage. This script only has to split those residuals by where the insertion's carriers come from.

GROUPS, from the founder panel (231 founders, origin-bio1 tertiles of 77):
  cold-specific  >= 2 cold carriers and ZERO warm carriers
  warm-specific  >= 2 warm carriers and ZERO cold carriers
  shared         carriers in both

Env: kmate. Read-only apart from a printed report.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
STORE = lib.AF_STORE
MIN_P0, SV_BP, MIN_MAC = 0.02, 50, 12


def main():
    F = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    grp = F.group.to_numpy()
    idx = np.load(f"{STORE}/index_nonsnp.npz")
    ch = idx["chrom"].astype("U5"); pos = idx["pos"].astype(np.int64)
    dl = idx["alt_len"].astype(np.int64) - idx["ref_len"].astype(np.int64)
    p0 = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
    keep = (lib.founder_panel_keep(ch, pos, min_mac=MIN_MAC) & (p0 >= MIN_P0) & (p0 <= 1 - MIN_P0))
    cols = np.where(keep)[0]
    is_ins = dl[cols] > SV_BP
    ins_cols = cols[is_ins]                                   # SAME order as res_ins / p0_ins
    print(f"{ins_cols.size:,} tested insertions")

    n = {g: np.zeros(ins_cols.size, np.int64) for g in ("cold", "temperate", "warm")}
    off = 0
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        prl = meta["ref_len"].astype(np.int64); pal = meta["alt_len"].astype(np.int64)
        nonsnp = np.where(~((prl == 1) & (pal == 1)))[0]
        nrow = int((ch == C).sum())
        sel = (ins_cols >= off) & (ins_cols < off + nrow)
        if sel.any():
            vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
            P = (vp[:, nonsnp[ins_cols[sel] - off]] > 0).astype(np.float64)
            for g in n:
                n[g][sel] = np.asarray(P[grp == g].sum(0)).ravel().astype(np.int64)
        off += nrow
        print(f"  [{C}] done", flush=True)

    cold_sp = (n["cold"] >= 2) & (n["warm"] == 0)
    warm_sp = (n["warm"] >= 2) & (n["cold"] == 0)
    shared = (n["cold"] > 0) & (n["warm"] > 0)
    print(f"\ncold-specific {cold_sp.sum():,} | warm-specific {warm_sp.sum():,} | shared {shared.sum():,}"
          f" | other {int(ins_cols.size - cold_sp.sum() - warm_sp.sum() - shared.sum()):,}")

    P = np.load(f"{OUT}/sv_temporal_climate_matched_persite.npz")
    D = pd.read_csv(f"{OUT}/sv_temporal_climate_matched_sites.csv").sort_values("bio1")
    GRPS = (("cold-specific", cold_sp), ("shared", shared), ("warm-specific", warm_sp))

    rows = []
    for _, r in D.iterrows():
        res = P[f"{int(r.site)}_res_ins_matched"]
        d = dict(site=int(r.site), bio1=r.bio1, bio18=r.bio18)
        for lab, m in GRPS:
            v = res[m]; v = v[~np.isnan(v)]
            d[lab] = float(np.median(v)) if v.size >= 20 else np.nan
        rows.append(d)
    R = pd.DataFrame(rows)
    R.to_csv(f"{OUT}/sv_cold_specific_global.csv", index=False)

    print("\nmedian residual per garden (vs SNPs matched on p0 + carrier climate), by carrier origin")
    print("if insertions are simply deleterious, EVERY row should be negative in EVERY garden group\n")
    R["terc"] = pd.qcut(R.bio1, 3, labels=["cold gardens", "mid gardens", "hot gardens"])
    print(f"{'insertion group':<16}{'garden group':<15}{'median':>10}{'<0 at':>9}{'wilcoxon p':>12}")
    for lab, _ in GRPS:
        for g in ("cold gardens", "mid gardens", "hot gardens"):
            v = R[R.terc == g][lab].dropna()
            if v.size < 3:
                continue
            print(f"{lab:<16}{g:<15}{v.median():>+10.4f}{f'{int((v<0).sum())}/{len(v)}':>9}"
                  f"{stats.wilcoxon(v).pvalue:>12.3f}")
        v = R[lab].dropna()
        print(f"{lab:<16}{'ALL GARDENS':<15}{v.median():>+10.4f}{f'{int((v<0).sum())}/{len(v)}':>9}"
              f"{stats.wilcoxon(v).pvalue:>12.3f}")
        r1 = stats.spearmanr(R.dropna(subset=[lab]).bio1, v)
        print(f"{'':<16}{'vs bio1':<15}{'rho':>10} {r1.statistic:+.3f}   p={r1.pvalue:.4f}\n")
    print(f"[wrote] {OUT}/sv_cold_specific_global.csv")


if __name__ == "__main__":
    main()
