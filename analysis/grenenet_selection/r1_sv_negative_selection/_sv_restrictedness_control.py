#!/usr/bin/env python
"""Is the global decline a property of INSERTIONS, or of being LINEAGE-RESTRICTED?

Cold-specific insertions decline in 29/31 gardens -- but so do warm-specific ones (21/31), with no
climate gradient, while shared insertions do not decline at all. The property that predicts decline
is therefore restriction to one origin tertile, not cold-carriage. Two readings:

  (a) lineage-restricted insertions really are younger / more deleterious; or
  (b) restrictedness itself manufactures a decline -- a variant confined to one clade has its pool
      frequency read off a narrower founder projection, and systematic bias there becomes a
      consistently negative logit slope.

The control: apply the IDENTICAL cold-specific / warm-specific / shared split to INDELS, which are
the reference class, and ask whether restricted indels also decline globally. If they do, (b) holds
and the cold-specific insertion result dissolves. If they do not, the effect is specific to
insertions.

Reproduces the indel subsample of `_sv_temporal_climate_matched.py` exactly (same seed, same steps).

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
MIN_P0, SV_BP, MIN_MAC, NSUB_IND = 0.02, 50, 12, 25000


def main():
    F = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    grp = F.group.to_numpy()
    idx = np.load(f"{STORE}/index_nonsnp.npz")
    ch = idx["chrom"].astype("U5"); pos = idx["pos"].astype(np.int64)
    dl = idx["alt_len"].astype(np.int64) - idx["ref_len"].astype(np.int64)
    p0 = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
    keep = (lib.founder_panel_keep(ch, pos, min_mac=MIN_MAC) & (p0 >= MIN_P0) & (p0 <= 1 - MIN_P0))
    cols = np.where(keep)[0]
    dln = dl[cols]; is_ins = dln > SV_BP; is_ind = (np.abs(dln) >= 1) & (np.abs(dln) <= SV_BP)

    # carrier counts per origin tertile, for every tested non-SNP record
    n = {g: np.zeros(cols.size, np.int64) for g in ("cold", "temperate", "warm")}
    tot = np.zeros(cols.size, np.int64)
    off = 0
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        prl = meta["ref_len"].astype(np.int64); pal = meta["alt_len"].astype(np.int64)
        nonsnp = np.where(~((prl == 1) & (pal == 1)))[0]
        nrow = int((ch == C).sum())
        sel = (cols >= off) & (cols < off + nrow)
        if sel.any():
            vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
            P = (vp[:, nonsnp[cols[sel] - off]] > 0).astype(np.float64)
            for g in n:
                n[g][sel] = np.asarray(P[grp == g].sum(0)).ravel().astype(np.int64)
            tot[sel] = np.asarray(P.sum(0)).ravel().astype(np.int64)
        off += nrow
        print(f"  [{C}] done", flush=True)

    ok = tot > 0                                     # == ~isnan(carrier climate) in the source script
    rng = np.random.default_rng(0)
    ii = np.where(is_ind & ok)[0]
    i_sub = np.sort(rng.choice(ii, min(NSUB_IND, ii.size), replace=False))
    print(f"\nreproduced indel subsample: {i_sub.size:,} of {ii.size:,}")

    def split(mask_idx):
        c, w = n["cold"][mask_idx], n["warm"][mask_idx]
        return (("cold-specific", (c >= 2) & (w == 0)),
                ("shared", (c > 0) & (w > 0)),
                ("warm-specific", (w >= 2) & (c == 0)))

    P = np.load(f"{OUT}/sv_temporal_climate_matched_persite.npz")
    D = pd.read_csv(f"{OUT}/sv_temporal_climate_matched_sites.csv").sort_values("bio1")
    sites = D.site.astype(int).tolist()

    print(f"\n{'class':<10}{'group':<16}{'n vars':>8}{'median across gardens':>23}{'<0 at':>9}{'wilcoxon p':>12}{'rho bio1':>11}")
    is_del = -dln > SV_BP
    MED = {}
    for cls, key, midx in (("INDEL", "indel", i_sub), ("INS", "ins", np.where(is_ins)[0]),
                           ("DEL", "del", np.where(is_del)[0])):
        for lab, m in split(midx):
            if m.sum() < 30:
                print(f"{cls:<10}{lab:<16}{int(m.sum()):>8}  (too few)")
                continue
            med = []
            for s in sites:
                v = P[f"{s}_res_{key}_matched"][m]; v = v[~np.isnan(v)]
                med.append(np.median(v) if v.size >= 20 else np.nan)
            med = np.array(med); g = ~np.isnan(med)
            MED[(cls, lab)] = med
            r = stats.spearmanr(D.bio1.to_numpy()[g], med[g])
            print(f"{cls:<10}{lab:<16}{int(m.sum()):>8}{np.median(med[g]):>+23.4f}"
                  f"{f'{int((med[g]<0).sum())}/{g.sum()}':>9}{stats.wilcoxon(med[g]).pvalue:>12.4f}"
                  f"{r.statistic:>+11.3f}")
        print()

    # is the cold-specific decline deeper than the warm-specific one? paired across gardens
    print("paired across gardens, cold-specific vs warm-specific (same gardens, same matching):")
    for cls in ("INS", "DEL", "INDEL"):
        a = MED.get((cls, "cold-specific")); b = MED.get((cls, "warm-specific"))
        if a is None or b is None:
            continue
        g = ~np.isnan(a) & ~np.isnan(b)
        d = a[g] - b[g]
        print(f"  {cls:<6} cold-warm = {np.median(d):+.4f}  cold deeper at {int((d<0).sum())}/{g.sum()}"
              f"  wilcoxon p={stats.wilcoxon(d).pvalue:.4f}")


if __name__ == "__main__":
    main()
