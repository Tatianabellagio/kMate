#!/usr/bin/env python
"""Does the temporal SV depletion survive matching SNPs on ECOTYPE-OF-ORIGIN CLIMATE as well as p0?

THE CONFOUND. The per-site SV-vs-SNP shift is de-trended against SNPs matched on starting frequency
p0. But pool AF is a founder projection, so when a garden's ecotype composition shifts (warm-origin
ecotypes winning at hot gardens), every variant carried preferentially by cold-origin ecotypes falls.
Insertion alleles are carried slightly more often by cold-origin ecotypes (per-founder count
rho = -0.36 with origin bio1), so they fall relative to an unmatched SNP average -- with no selection
on insertions at all. Frequency matching cannot remove this: two variants at the same p0 can sit on
completely different ecotype backgrounds, and it is the background that garden turnover acts on.

THE FIX. Match the SNP comparison set on BOTH p0 and carrier-origin climate, so the contrast is taken
between variants carried by ecotypes of the same origin climate. carrier_bio1 = mean origin bio1 of
the founders carrying the variant, computed from the founder panel for SVs and SNPs alike.

READ THE CONTROL FIRST. The p0-only shift is recomputed here in the same run. If it does not land on
the published rho (~ -0.39 bio1 / +0.73 bio18 for insertions), the matched number means nothing.

SV_UNIT=event: SVs scored as Truvari EVENTS (`_build_event_units.py`): pool AF of an event = sum of its
alleles' AFs (capped at 1), p0 likewise; filtered on the event's own union carriers (MIN_MAC..N-MIN_MAC)
and call rate (>= 90%); carrier climate = mean bio1 of the union carriers. Short indels stay per record.
Outputs tagged `_event{MIN_MAC}`.
Env: kmate. Writes results/sv_adaptive/sv_temporal_climate_matched_{sites,summary}.csv.
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
SV_BP, NBIN1, NBINP, NBINC, MINCELL = 50, 10, 10, 5, 20
# SV_MIN_P0: seed-mix frequency floor, applied symmetrically (p0 >= MIN_P0 and <= 1 - MIN_P0) to the
# tested variants AND the SNP comparison set. 0.02 = the published analysis. Lowering it admits
# low-frequency variants on both sides; the p0 bins are quantiles of the combined set, so they adapt.
# Promoted to an env var 2026-09-24 (was hardcoded); default unchanged, so existing runs reproduce.
MIN_P0 = float(os.environ.get("SV_MIN_P0", 0.02))
MIN_MAC = int(os.environ.get("SV_MIN_MAC", 12))     # founder carrier floor; 12 = the published analysis
TAG = ("" if MIN_MAC == 12 else f"_mac{MIN_MAC}")      # output suffix so runs do not overwrite each other
CALL_MIN = float(os.environ.get("SV_CALL_MIN", 0.9))   # founder call-rate floor; != 0.9 -> outputs tagged _cr{pct}
CRTAG = "" if CALL_MIN == 0.9 else f"_cr{int(round(CALL_MIN * 100))}"
FILTER = os.environ.get("SV_FILTER", "position")   # "record" = lib.founder_record_keep (each allele on its OWN carriers)
def panel_keep(kind, chrom, pos, min_mac):
    """carrier/call-rate floor. 'position' reproduces the published analysis (lib.founder_panel_keep:
    rare alleles pass on a common sibling at the same position); 'record' judges each allele alone."""
    if FILTER == "record":
        return lib.founder_record_keep(kind, np.arange(len(pos)), min_mac=min_mac, called_min=CALL_MIN)
    return lib.founder_panel_keep(chrom, pos, min_mac=min_mac)
if FILTER == "record": TAG = f"_rec{MIN_MAC}"
UNIT = os.environ.get("SV_UNIT", "allele")         # "event": SVs as Truvari events (`_build_event_units.py`);
if UNIT == "event":                                # short indels stay per record; filter on the unit's OWN
    FILTER = "record"; TAG = f"_event{MIN_MAC}"    # union carriers / call rate; pool AF = sum of member AFs
if CRTAG:
    assert FILTER == "record", "SV_CALL_MIN only applies to the record / event filter"
    TAG += CRTAG
P0TAG = "" if abs(MIN_P0 - 0.02) < 1e-12 else f"_p0{MIN_P0:g}".replace(".", "")
TAG += P0TAG
NSUB_IND = 25000            # indels stored subsampled, for the per-site grid
SNP_SUB = int(os.environ.get("SV_SNP_SUB", 0))   # >0: also store residuals for a random SNP subsample
                                                 # (per-context SNP reference; own rng, so the indel draw is unchanged)


def carrier_climate(kind, cols):
    """mean origin bio1 of the founders carrying each record, in af_store column order."""
    F = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    b1 = F.bio1.to_numpy(np.float64)
    idx = np.load(f"{STORE}/index_{kind}.npz"); ch = idx["chrom"].astype("U5")
    rl = idx["ref_len"].astype(np.int64); al = idx["alt_len"].astype(np.int64)
    out = np.full(ch.size, np.nan); off = 0
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        prl = meta["ref_len"].astype(np.int64); pal = meta["alt_len"].astype(np.int64)
        issnp = (prl == 1) & (pal == 1)
        sel = np.where(issnp if kind == "snp" else ~issnp)[0]
        m = ch == C; n = int(m.sum())
        assert n == sel.size, f"{C}/{kind}: af_store {n:,} vs panel {sel.size:,}"
        assert np.array_equal(prl[sel], rl[m]) and np.array_equal(pal[sel], al[m]), f"{C}/{kind}: order"
        loc = cols[(cols >= off) & (cols < off + n)] - off           # wanted, chrom-local
        if loc.size:
            vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
            P = (vp[:, sel[loc]] > 0).astype(np.float64)             # [231 x k]
            nc = np.asarray(P.sum(0)).ravel()
            out[off + loc] = np.asarray((b1[None, :] @ P)).ravel() / np.maximum(nc, 1)
        off += n
        print(f"    [{C}/{kind}] {loc.size:,} records", flush=True)
    return out


def cells(p0_t, c_t, p0_n, c_n, joint):
    """Bin target (non-SNP) and SNP records the same way. joint=False: p0 bins only (the published
    de-trending). joint=True: p0 x carrier-climate bins, so the SNPs a given record is compared
    against are carried by ecotypes of the same origin climate. Both use the SAME p0 resolution,
    so the only difference between them is the climate dimension."""
    ep = np.quantile(np.r_[p0_t, p0_n], np.linspace(0, 1, NBINP + 1))[1:-1]
    if not joint:
        return np.digitize(p0_t, ep), np.digitize(p0_n, ep), NBINP
    ec = np.quantile(np.r_[c_t, c_n], np.linspace(0, 1, NBINC + 1))[1:-1]
    return (np.digitize(p0_t, ep) * NBINC + np.digitize(c_t, ec),
            np.digitize(p0_n, ep) * NBINC + np.digitize(c_n, ec), NBINP * NBINC)


def residual(s_t, bt, s_n, bn, nb):
    """each record's s minus the median s of the SNPs in its OWN cell (NaN if the cell is thin)."""
    med = np.full(nb, np.nan)
    for b in range(nb):
        c = s_n[bn == b]
        if c.size >= MINCELL:
            med[b] = np.median(c)
    return s_t - med[bt]


def shift_from(res, bt, nb):
    """mean over cells of the median residual -- same statistic as before, built from residuals."""
    d = [np.median(res[(bt == b) & ~np.isnan(res)]) for b in range(nb)
         if int(((bt == b) & ~np.isnan(res)).sum()) >= MINCELL]
    return float(np.mean(d)) if d else np.nan


def plot_replicate_sz_units(site, A, need, p0c):
    """as _temporal_s_plots_snp_vs_nonsnp.plot_replicate_sz, for UNITS: each unit's pool AF is
    A @ af[need] (sum of its member records' AFs), capped at 1."""
    L0 = np.log(np.clip(p0c, pmod_EPS, 1 - pmod_EPS) / (1 - np.clip(p0c, pmod_EPS, 1 - pmod_EPS)))
    plot_pts = {}
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
            af = np.minimum(A @ np.asarray(mat[ridx])[need].astype(np.float64), 1.0)
            plot_pts.setdefault(int(plot), []).append(
                (float(g), np.log(np.clip(af, pmod_EPS, 1 - pmod_EPS) / (1 - np.clip(af, pmod_EPS, 1 - pmod_EPS)))))
    slopes = []
    for plot, pts in plot_pts.items():
        ts = np.array([0.0] + [t_ for t_, _ in pts])
        if ts.size < 2:
            continue
        Y = np.vstack([L0] + [y for _, y in pts])
        tc = ts - ts.mean(); slopes.append((tc @ Y) / float((tc ** 2).sum()))
    S = np.vstack(slopes); n = S.shape[0]
    return S.mean(0).astype(np.float32), None, None, n


def main():
    import importlib.util
    _sp = importlib.util.spec_from_file_location(
        "pm", os.path.join(os.path.dirname(os.path.abspath(__file__)), "_temporal_s_plots_snp_vs_nonsnp.py"))
    pmod = importlib.util.module_from_spec(_sp); _sp.loader.exec_module(pmod)
    global pmod_EPS; pmod_EPS = pmod.EPS

    idx_non = np.load(f"{STORE}/index_nonsnp.npz"); idx_snp = np.load(f"{STORE}/index_snp.npz")
    ch_non = idx_non["chrom"].astype("U5"); pos_non = idx_non["pos"].astype(np.int64)
    dl = idx_non["alt_len"].astype(np.int64) - idx_non["ref_len"].astype(np.int64)
    p0_non = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
    p0_snp = np.load(f"{STORE}/p0_snp.npy").astype(np.float64)

    keep_non = (panel_keep("nonsnp", ch_non, pos_non, MIN_MAC)
                & (p0_non >= MIN_P0) & (p0_non <= 1 - MIN_P0))
    keep_snp = (panel_keep("snp", idx_snp["chrom"].astype("U5"), idx_snp["pos"].astype(np.int64), MIN_MAC)
                & (p0_snp >= MIN_P0) & (p0_snp <= 1 - MIN_P0))
    cols_snp = np.where(keep_snp)[0]
    if UNIT == "event":
        EU = np.load(f"{OUT}/nonsnp_event_units.npz"); NF = int(EU["n_founders"])
        A_all = sp.csr_matrix((np.ones(EU["indices"].size), EU["indices"], EU["indptr"]),
                              shape=(EU["rep"].size, dl.size))                  # units x af_store rows
        p0_u = np.minimum(A_all @ p0_non, 1.0)
        keep_u = ((EU["n"] >= MIN_MAC) & (EU["n"] <= NF - MIN_MAC) & (EU["ncal"] >= CALL_MIN * NF)
                  & (p0_u >= MIN_P0) & (p0_u <= 1 - MIN_P0))
        cols_non = np.where(keep_u)[0]                                          # UNIT ids
        A_keep = A_all[cols_non]; need = np.unique(A_keep.indices); A_small = A_keep[:, need].tocsr()
        dln = EU["dlen"][cols_non]; p0n = p0_u[cols_non]
    else:
        cols_non = np.where(keep_non)[0]
        dln = dl[cols_non]; p0n = p0_non[cols_non]
    is_ins = dln > SV_BP; is_del = -dln > SV_BP; is_sv = is_ins | is_del
    print(f"[1] {cols_non.size:,} non-SNP {'units' if UNIT == 'event' else 'records'} "
          f"({is_ins.sum():,} ins, {is_del.sum():,} del), {cols_snp.size:,} SNP")

    print("[2] carrier-origin climate")
    c_non = EU["b1"][cols_non] if UNIT == "event" else carrier_climate("nonsnp", cols_non)[cols_non]
    c_snp = carrier_climate("snp", cols_snp)[cols_snp]
    p0s = p0_snp[cols_snp]
    ok = ~np.isnan(c_non); oks = ~np.isnan(c_snp)
    print(f"    carrier bio1: SV {np.nanmean(c_non[is_sv]):.2f} C, SNP {np.nanmean(c_snp):.2f} C, "
          f"insertions {np.nanmean(c_non[is_ins]):.2f} C, deletions {np.nanmean(c_non[is_del]):.2f} C")
    print(f"    (if insertions sit colder than SNPs, composition turnover alone depletes them)")

    cmeta = pd.concat([pd.read_csv(f"{PM}/pool_gen{g}_snp.meta.csv") for g in (1, 2, 3)],
                      ignore_index=True).groupby("site")[["bio1", "bio18"]].mean()
    sites = sorted({int(s) for f in glob.glob(f"{PM}/pool_gen*_snp.meta.csv")
                    for s in pd.read_csv(f).site.unique()})

    is_ind = (np.abs(dln) >= 1) & (np.abs(dln) <= SV_BP)
    rng = np.random.default_rng(0)
    ii = np.where(is_ind & ok)[0]
    i_sub = np.sort(rng.choice(ii, min(NSUB_IND, ii.size), replace=False))   # indels, for the grid
    CLS = (("sv", is_sv), ("ins", is_ins), ("del", is_del))
    s_sub = (np.sort(np.random.default_rng(1).choice(np.where(oks)[0], min(SNP_SUB, int(oks.sum())), replace=False))
             if SNP_SUB > 0 else np.zeros(0, np.int64))              # positions within cols_snp

    print("[3] per site")
    rows = []; sub = {}
    for site in sites:
        s_sn, _, _, n = pmod.plot_replicate_sz(site, "snp", cols_snp, p0s)
        if n < 2 or site not in cmeta.index:
            continue
        s_no = (plot_replicate_sz_units(site, A_small, need, p0n)[0] if UNIT == "event"
                else pmod.plot_replicate_sz(site, "nonsnp", cols_non, p0n)[0])
        r = dict(site=site, n_plots=n, bio1=float(cmeta.loc[site, "bio1"]), bio18=float(cmeta.loc[site, "bio18"]))
        for tag, joint in (("p0only", False), ("matched", True)):
            # one binning for the whole non-SNP set, so SV/indel residuals share a baseline
            bt, bn, nb = cells(p0n[ok], c_non[ok], p0s[oks], c_snp[oks], joint)
            res_all = residual(s_no[ok], bt, s_sn[oks], bn, nb)
            for lab, msk in CLS:
                sel = msk[ok]
                r[f"{lab}_{tag}"] = shift_from(res_all[sel], bt[sel], nb)
            # per-variant residuals for the per-site grid
            full = np.full(cols_non.size, np.nan); full[ok] = res_all
            for lab, msk in CLS:
                sub[f"{site}_res_{lab}_{tag}"] = full[msk].astype(np.float32)
            sub[f"{site}_res_indel_{tag}"] = full[i_sub].astype(np.float32)
            if s_sub.size:
                # SNP residuals against the SAME cells: each SNP minus its own cell's SNP median
                ok_pos = np.where(oks)[0]; rank = np.searchsorted(ok_pos, s_sub)
                res_sn = residual(s_sn[oks], bn, s_sn[oks], bn, nb)
                sub[f"{site}_res_snp_{tag}"] = res_sn[rank].astype(np.float32)
        rows.append(r)
        print(f"  site {site:>2} bio1={r['bio1']:>5.1f}: ins p0only {r['ins_p0only']:+.4f} "
              f"-> matched {r['ins_matched']:+.4f}", flush=True)

    for lab, msk in CLS:
        sub[f"p0_{lab}"] = p0n[msk].astype(np.float32)
    sub["p0_indel"] = p0n[i_sub].astype(np.float32)
    sub["cols_non"] = cols_non.astype(np.int64)        # af_store rows, or UNIT ids when SV_UNIT=event
    sub["i_sub"] = i_sub.astype(np.int64)              # indel subsample, positions within cols_non
    if s_sub.size:
        sub["snp_rows"] = cols_snp[s_sub].astype(np.int64)            # af_store index_snp rows
        sub["p0_snp_sub"] = p0s[s_sub].astype(np.float32)
    D = pd.DataFrame(rows); D.to_csv(f"{OUT}/sv_temporal_climate_matched_sites{TAG}.csv", index=False)
    np.savez_compressed(f"{OUT}/sv_temporal_climate_matched_persite{TAG}.npz", **sub)
    print(f"[wrote] {OUT}/sv_temporal_climate_matched_persite{TAG}.npz")
    out = []
    print(f"\n=== per-site shift vs garden climate ({len(D)} sites) ===")
    print(f"{'measure':<16}{'median shift':>14}{'<0 at':>9}{'rho bio1':>11}{'p':>9}{'rho bio18':>12}{'p':>9}")
    for lab in ("sv", "ins", "del"):
        for tag in ("p0only", "matched"):
            c = f"{lab}_{tag}"; d = D[["bio1", "bio18", c]].dropna()
            r1 = stats.spearmanr(d.bio1, d[c]); r2 = stats.spearmanr(d.bio18, d[c])
            print(f"{c:<16}{d[c].median():>+14.4f}{f'{(d[c]<0).sum()}/{len(d)}':>9}"
                  f"{r1.statistic:>+11.3f}{r1.pvalue:>9.4f}{r2.statistic:>+12.3f}{r2.pvalue:>9.4f}")
            out.append(dict(measure=c, n_sites=len(d), median_shift=d[c].median(),
                            n_neg=int((d[c] < 0).sum()), rho_bio1=r1.statistic, p_bio1=r1.pvalue,
                            rho_bio18=r2.statistic, p_bio18=r2.pvalue))
    pd.DataFrame(out).to_csv(f"{OUT}/sv_temporal_climate_matched_summary{TAG}.csv", index=False)
    print(f"\n[wrote] {OUT}/sv_temporal_climate_matched_{{sites,summary}}{TAG}.csv")


if __name__ == "__main__":
    main()
