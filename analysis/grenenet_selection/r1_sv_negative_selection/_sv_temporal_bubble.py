#!/usr/bin/env python
"""Does the temporal SV-purging climate gradient survive counting BUBBLES instead of ALLELES?

`temporal_s_consolidated.ipynb` reports (i) a net de-trended shift of ~0 and (ii) a climate gradient
in the per-site shift: rho ~ -0.48 with bio1, ~ +0.65 with bio18, carried by INSERTIONS. Those are
computed on biallelic panel records = ALLELES. arch3 is the biallelic representation of a
multi-allelic graph, so one insertion event appears once per distinct inserted sequence; 83% of the
tested insertions sit in a bubble with other alleles and 38% share a bubble with another TESTED
allele, so the tested variants are not independent and can tag lineages rather than events.

s is a logit slope, so it cannot be re-aggregated from per-allele s values (the slope of a sum is
not the sum of slopes). It has to be recomputed from the frequency trajectories:

    bubble AF in a pool = SUM of its alleles' AFs   (alleles of one bubble are mutually exclusive
                                                     on a haplotype, so they add)
    p0(bubble)          = SUM of its alleles' p0
    then the SAME plot-replicate logit-slope estimator as `_temporal_s_plots_snp_vs_nonsnp.py`.

Aggregate FIRST, logit SECOND. Both units are computed in one pass over the pool matrices so they
see identical data, and the ALLELE result is reproduced here as a control: if it does not land on
the published rho the bubble comparison means nothing.

af_store's non-SNP index is row-for-row identical to the arch3 panel's non-SNP records on all five
chromosomes (verified), so record -> bubble is a positional lookup, not a key join -- the
multiallelic chrom:pos trap does not apply.

Env: kmate. Writes results/sv_adaptive/sv_temporal_bubble_{sites,summary}.csv.
"""
import os, sys, glob
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

os.chdir("/global/scratch/users/tbellg/kmate")
STORE = lib.AF_STORE
PM = f"{lib.GEA}/common/results/pool_matrices"
OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
MIN_P0, SV_BP, MIN_MAC, CALL_MIN, NBIN, EPS = 0.02, 50, 12, 0.9, 10, 1e-3


def logit(p):
    q = np.clip(p, EPS, 1 - EPS)
    return np.log(q / (1 - q))


def build_bubbles():
    """SV allele columns of af_store, their bubble grouping, class and founder-panel stats."""
    idx = np.load(f"{STORE}/index_nonsnp.npz")
    ch = idx["chrom"].astype("U5")
    rl = idx["ref_len"].astype(np.int64); al = idx["alt_len"].astype(np.int64); dl = al - rl
    BM = np.load(f"{OUT}/arch3_bubble_map.npz")

    cols, bub, ins = [], [], []
    n_alt, n_cal = [], []
    off = 0
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        prl = meta["ref_len"].astype(np.int64); pal = meta["alt_len"].astype(np.int64)
        nonsnp = np.where(~((prl == 1) & (pal == 1)))[0]          # panel index of each af_store row
        m = ch == C
        assert int(m.sum()) == nonsnp.size, f"{C}: af_store/panel non-SNP count differs"
        assert np.array_equal(prl[nonsnp], rl[m]) and np.array_equal(pal[nonsnp], al[m]), f"{C}: order"
        rec2bub = pd.Series(BM[f"{C}_bubble"], index=BM[f"{C}_rec"])
        loc = np.where(np.abs(dl[m]) > SV_BP)[0]                  # SV rows of this chrom, af_store-local
        b = rec2bub.reindex(nonsnp[loc]).to_numpy()
        assert not np.isnan(b).any() and (b >= 0).all(), f"{C}: SV rows without a bubble code"
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        pc = nonsnp[loc]
        n_alt.append(np.asarray(vp[:, pc].sum(0)).ravel()); n_cal.append(np.asarray(vc[:, pc].sum(0)).ravel())
        cols.append(off + loc); bub.append(b.astype(np.int64) + 10**7 * lib.CHROMS.index(C))
        ins.append(dl[m][loc] > 0)
        off += int(m.sum())
        print(f"  [{C}] {loc.size:,} SV alleles", flush=True)

    cols = np.concatenate(cols); bub = np.concatenate(bub); ins = np.concatenate(ins)
    n_alt = np.concatenate(n_alt); n_cal = np.concatenate(n_cal)
    ub, li = np.unique(bub, return_inverse=True); li = li.ravel()
    M = sp.csr_matrix((np.ones(cols.size), (li, np.arange(cols.size))), shape=(ub.size, cols.size))
    k = np.bincount(li, minlength=ub.size)                        # alleles per bubble
    nins = np.bincount(li, weights=ins.astype(float), minlength=ub.size)
    bcls = np.where(nins == k, "ins", np.where(nins == 0, "del", "mixed"))
    return cols, M, li, k, bcls, ins, n_alt, n_cal


def exclusive_mask(cols, li, nbub):
    """Is a bubble's allele set mutually EXCLUSIVE on a haplotype?

    Summing allele frequencies to get a bubble frequency is only valid if no haplotype carries two
    of the bubble's alleles. arch3 records are ATOMIC variants nested in a bubble, so a single path
    through the bubble can carry several of them -- and it does: summed p0 exceeds 1 for 2.6% of
    bubbles (max 17.3), concentrated in high-allele bubbles. Only 55.9% of k >= 10 bubbles are
    exclusive, against 100% of k = 1. Tested here directly on the founder panel: a bubble is
    exclusive if no founder carries more than one of its records. Non-exclusive bubbles must be
    dropped, not clipped -- their summed frequency is not a frequency."""
    idx = np.load(f"{STORE}/index_nonsnp.npz"); ch = idx["chrom"].astype("U5")
    mx = np.zeros(nbub)
    off = 0
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        prl = meta["ref_len"].astype(np.int64); pal = meta["alt_len"].astype(np.int64)
        nonsnp = np.where(~((prl == 1) & (pal == 1)))[0]
        nrow = int((ch == C).sum())
        sel = (cols >= off) & (cols < off + nrow)
        pc = nonsnp[cols[sel] - off]
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
        P = (vp[:, pc] > 0).astype(np.float64)
        L = sp.csr_matrix((np.ones(int(sel.sum())), (li[sel], np.arange(int(sel.sum())))),
                          shape=(nbub, int(sel.sum())))
        cnt = np.asarray((L @ P.T.tocsr()).todense())      # [nbub x 231] records carried per founder
        mx = np.maximum(mx, cnt.max(1))
        off += nrow
    return mx <= 1


def bubble_founder_stats(cols, li, nbub):
    """carriers/called per bubble = union over its alleles, across the 231 founders."""
    idx = np.load(f"{STORE}/index_nonsnp.npz"); ch = idx["chrom"].astype("U5")
    rl = idx["ref_len"].astype(np.int64); al = idx["alt_len"].astype(np.int64)
    c = np.zeros(nbub, np.int64); m = np.zeros(nbub, np.int64)
    off = 0
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        prl = meta["ref_len"].astype(np.int64); pal = meta["alt_len"].astype(np.int64)
        nonsnp = np.where(~((prl == 1) & (pal == 1)))[0]
        nrow = int((ch == C).sum())
        sel = (cols >= off) & (cols < off + nrow)
        pc = nonsnp[cols[sel] - off]
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        P = (vp[:, pc] > 0).astype(np.int8).T.tocsr(); Cc = (vc[:, pc] > 0).astype(np.int8).T.tocsr()
        L = sp.csr_matrix((np.ones(sel.sum()), (li[sel], np.arange(sel.sum()))), shape=(nbub, int(sel.sum())))
        c += np.asarray(((L @ P) > 0).sum(1)).ravel().astype(np.int64)
        m += np.asarray(((L @ Cc) > 0).sum(1)).ravel().astype(np.int64)
        off += nrow
    return c, m


def site_slopes(site, cols, M, p0a, p0b, cols_ind=None, p0i=None):
    """One pass over a site's pools: per-plot logit slopes for SV ALLELES, SV BUBBLES and indels.

    Indels are read from the same pool rows so the three classes see identical data. They have no
    bubble aggregation available (the bubble map only covers records with |dlen| > 50 bp), so they
    stay at the allele unit and serve as the reference class in the per-site grid."""
    L0a, L0b = logit(p0a), logit(p0b)
    L0i = logit(p0i) if cols_ind is not None else None
    pts = {}
    for g in (1, 2, 3):
        try:
            meta = pd.read_csv(f"{PM}/pool_gen{g}_nonsnp.meta.csv")
        except FileNotFoundError:
            continue
        m4 = meta[meta.site == site]
        if len(m4) == 0:
            continue
        mat = np.load(f"{PM}/pool_gen{g}_nonsnp_af.npy", mmap_mode="r")
        for ridx, plot in zip(m4.index.to_numpy(), m4["plot"].to_numpy()):
            row = np.asarray(mat[ridx])
            af = row[cols].astype(np.float64)
            li_ = logit(row[cols_ind].astype(np.float64)) if cols_ind is not None else None
            pts.setdefault(int(plot), []).append((float(g), logit(af), logit(M @ af), li_))
    SA, SB, SI = [], [], []
    for plot, p in pts.items():
        ts = np.array([0.0] + [t for t, _, _, _ in p])
        if ts.size < 2:
            continue
        tc = ts - ts.mean(); sst = float((tc ** 2).sum())
        SA.append((tc @ np.vstack([L0a] + [a for _, a, _, _ in p])) / sst)
        SB.append((tc @ np.vstack([L0b] + [b for _, _, b, _ in p])) / sst)
        if cols_ind is not None:
            SI.append((tc @ np.vstack([L0i] + [c for _, _, _, c in p])) / sst)
    if not SA:
        return None, None, None, 0
    si = np.vstack(SI).mean(0) if SI else None
    return np.vstack(SA).mean(0), np.vstack(SB).mean(0), si, len(SA)


def main():
    import importlib.util
    _sp = importlib.util.spec_from_file_location(
        "pm", os.path.join(os.path.dirname(os.path.abspath(__file__)), "_temporal_s_plots_snp_vs_nonsnp.py"))
    pmod = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(pmod)

    print("[1] bubbles")
    cols, M, li, k, bcls, ins, n_alt, n_cal = build_bubbles()
    nbub = M.shape[0]
    p0_non = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
    p0a = p0_non[cols]; p0b = np.asarray(M @ p0a).ravel()
    cb, mb = bubble_founder_stats(cols, li, nbub)
    excl = exclusive_mask(cols, li, nbub)
    print(f"    exclusive bubbles (summation valid): {excl.sum():,} of {nbub:,} ({100*excl.mean():.1f}%)")
    print(f"    {cols.size:,} SV alleles -> {nbub:,} bubbles "
          f"({(bcls=='ins').sum():,} ins, {(bcls=='del').sum():,} del, {(bcls=='mixed').sum():,} mixed)")

    # masks. alleles: the existing per-record floor. bubbles: the same floor applied to the EVENT.
    idx = np.load(f"{STORE}/index_nonsnp.npz")
    keep_a = (lib.founder_panel_keep(idx["chrom"].astype("U5")[cols], idx["pos"].astype(np.int64)[cols],
                                     min_mac=MIN_MAC) & (p0a >= MIN_P0) & (p0a <= 1 - MIN_P0))
    keep_b = ((cb >= MIN_MAC) & (cb <= 231 - MIN_MAC) & (mb >= CALL_MIN * 231)
              & (p0b >= MIN_P0) & (p0b <= 1 - MIN_P0) & excl)   # summation is only valid if exclusive
    print(f"    tested: {keep_a.sum():,} alleles, {keep_b.sum():,} bubbles "
          f"({(keep_b & (bcls=='ins')).sum():,} ins, {(keep_b & (bcls=='del')).sum():,} del)")
    print(f"    of tested alleles, share in a multi-allele bubble: {100*(k[li][keep_a]>1).mean():.1f}%")

    idx_snp = np.load(f"{STORE}/index_snp.npz")
    p0_snp = np.load(f"{STORE}/p0_snp.npy").astype(np.float64)
    ksnp = (lib.founder_panel_keep(idx_snp["chrom"].astype("U5"), idx_snp["pos"].astype(np.int64),
                                   min_mac=MIN_MAC) & (p0_snp >= MIN_P0) & (p0_snp <= 1 - MIN_P0))
    cols_snp = np.where(ksnp)[0]; p0s = p0_snp[cols_snp]
    p0q = np.quantile(p0s, np.linspace(0, 1, NBIN + 1)); p0q[-1] += 1e-6
    binf = lambda p: np.clip(np.digitize(p, p0q[1:-1]), 0, NBIN - 1)
    bsnp = binf(p0s); ba = binf(p0a); bb = binf(p0b)

    # indels (1-50 bp) as the reference class for the per-site grid; subsampled to 25k like the
    # source analysis, with the same seed, so the grid's indel line is comparable to the published one
    rng = np.random.default_rng(0)
    dl_non = (np.load(f"{STORE}/index_nonsnp.npz")["alt_len"].astype(np.int64)
              - np.load(f"{STORE}/index_nonsnp.npz")["ref_len"].astype(np.int64))
    adl = np.abs(dl_non)
    ind_ok = np.where((adl >= 1) & (adl <= SV_BP)
                      & lib.founder_panel_keep(idx["chrom"].astype("U5"), idx["pos"].astype(np.int64),
                                               min_mac=MIN_MAC)
                      & (p0_non >= MIN_P0) & (p0_non <= 1 - MIN_P0))[0]
    cols_ind = np.sort(rng.choice(ind_ok, min(25000, ind_ok.size), replace=False))
    p0i = p0_non[cols_ind]
    print(f"    indel reference: {ind_ok.size:,} pass the filters, {cols_ind.size:,} sampled")
    i_sn = np.sort(rng.choice(cols_snp.size, min(25000, cols_snp.size), replace=False))

    clim = lib.load_climate()
    cmeta = pd.concat([pd.read_csv(f"{PM}/pool_gen{g}_snp.meta.csv") for g in (1, 2, 3)],
                      ignore_index=True).groupby("site")[["bio1", "bio18"]].mean()
    sites = sorted({int(s) for f in glob.glob(f"{PM}/pool_gen*_snp.meta.csv")
                    for s in pd.read_csv(f).site.unique()})

    print("[2] per site")
    rows = []; sub = {}
    for site in sites:
        s_sn, _, _, n = pmod.plot_replicate_sz(site, "snp", cols_snp, p0s)
        if n < 2 or site not in cmeta.index:
            continue
        sa, sb, si, n2 = site_slopes(site, cols, M, p0a, p0b, cols_ind, p0i)
        if sa is None:
            continue
        # per-variant (s, p0) for the per-site grid, one entry per class per site.
        # SV = BUBBLES here; SVallele keeps the allele unit for the side-by-side grid.
        sub[f"{site}_s_SNP"] = s_sn[i_sn].astype(np.float32); sub[f"{site}_p0_SNP"] = p0s[i_sn].astype(np.float32)
        sub[f"{site}_s_indel"] = si.astype(np.float32); sub[f"{site}_p0_indel"] = p0i.astype(np.float32)
        sub[f"{site}_s_SV"] = sb[keep_b].astype(np.float32); sub[f"{site}_p0_SV"] = p0b[keep_b].astype(np.float32)
        sub[f"{site}_s_SVallele"] = sa[keep_a].astype(np.float32)
        sub[f"{site}_p0_SVallele"] = p0a[keep_a].astype(np.float32)
        base = np.array([np.median(s_sn[bsnp == i]) if (bsnp == i).any() else np.nan for i in range(NBIN)])

        def shift(s, bins, mask):
            d = [np.median(s[mask & (bins == i)]) - base[i]
                 for i in range(NBIN) if (mask & (bins == i)).sum() >= 20]
            return float(np.mean(d)) if d else np.nan

        r = dict(site=site, n_plots=n, bio1=float(cmeta.loc[site, "bio1"]), bio18=float(cmeta.loc[site, "bio18"]),
                 allele_sv=shift(sa, ba, keep_a), bubble_sv=shift(sb, bb, keep_b),
                 allele_ins=shift(sa, ba, keep_a & ins), bubble_ins=shift(sb, bb, keep_b & (bcls == "ins")),
                 allele_del=shift(sa, ba, keep_a & ~ins), bubble_del=shift(sb, bb, keep_b & (bcls == "del")),
                 allele_ins_single=shift(sa, ba, keep_a & ins & (k[li] == 1)),
                 allele_ins_multi=shift(sa, ba, keep_a & ins & (k[li] > 1)),
                 # k=1: the allele IS the bubble, so these two must agree; a gap means the filters
                 # differ, not the units (the allele filter is position-based via founder_panel_keep)
                 allele_ins_k1=shift(sa, ba, keep_a & ins & (k[li] == 1)),
                 bubble_ins_k1=shift(sb, bb, keep_b & (bcls == "ins") & (k == 1)),
                 # the crux: multi-allele insertion EVENTS. their alleles carry the gradient; does
                 # the summed event carry it, or does it cancel within the bubble?
                 bubble_ins_multi=shift(sb, bb, keep_b & (bcls == "ins") & (k > 1)),
                 bubble_del_k1=shift(sb, bb, keep_b & (bcls == "del") & (k == 1)),
                 bubble_del_multi=shift(sb, bb, keep_b & (bcls == "del") & (k > 1)),
                 allele_del_k1=shift(sa, ba, keep_a & ~ins & (k[li] == 1)),
                 allele_del_multi=shift(sa, ba, keep_a & ~ins & (k[li] > 1)))
        if site == sites[0]:
            one = (k == 1)
            j = np.where(one)[0]
            a_of_b = np.asarray(M[j].argmax(1)).ravel()
            d = np.abs(sa[a_of_b] - sb[j])
            print(f"    [check] k=1 bubbles: max |s_allele - s_bubble| = {d.max():.2e} "
                  f"over {j.size:,} events (must be ~0)", flush=True)
        rows.append(r)
        print(f"  site {site:>2} bio1={r['bio1']:>5.1f}: allele {r['allele_sv']:+.4f} "
              f"bubble {r['bubble_sv']:+.4f} | ins allele {r['allele_ins']:+.4f} bubble {r['bubble_ins']:+.4f}",
              flush=True)

    D = pd.DataFrame(rows)
    D.to_csv(f"{OUT}/sv_temporal_bubble_sites.csv", index=False)
    np.savez_compressed(f"{OUT}/sv_temporal_bubble_persite.npz", **sub)
    print(f"[wrote] {OUT}/sv_temporal_bubble_persite.npz ({len(sub)//8} sites)")
    out = []
    print(f"\n=== per-site de-trended shift vs climate ({len(D)} sites) ===")
    print(f"{'measure':<20}{'median shift':>14}{'<0 at':>9}{'rho bio1':>12}{'p':>9}{'rho bio18':>12}{'p':>9}")
    for c in ("allele_sv", "bubble_sv", "allele_ins", "bubble_ins", "allele_del", "bubble_del",
              "allele_ins_k1", "bubble_ins_k1", "allele_ins_multi", "bubble_ins_multi",
              "allele_del_k1", "bubble_del_k1", "allele_del_multi", "bubble_del_multi",
              "allele_ins_single"):
        d = D[["bio1", "bio18", c]].dropna()
        r1 = stats.spearmanr(d.bio1, d[c]); r2 = stats.spearmanr(d.bio18, d[c])
        print(f"{c:<20}{d[c].median():>+14.4f}{f'{(d[c]<0).sum()}/{len(d)}':>9}"
              f"{r1.statistic:>+12.3f}{r1.pvalue:>9.4f}{r2.statistic:>+12.3f}{r2.pvalue:>9.4f}")
        out.append(dict(measure=c, n_sites=len(d), median_shift=d[c].median(), n_neg=int((d[c] < 0).sum()),
                        rho_bio1=r1.statistic, p_bio1=r1.pvalue, rho_bio18=r2.statistic, p_bio18=r2.pvalue))
    pd.DataFrame(out).to_csv(f"{OUT}/sv_temporal_bubble_summary.csv", index=False)
    print(f"\n[wrote] {OUT}/sv_temporal_bubble_{{sites,summary}}.csv")


if __name__ == "__main__":
    main()
