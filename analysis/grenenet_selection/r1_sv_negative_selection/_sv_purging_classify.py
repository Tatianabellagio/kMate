#!/usr/bin/env python
"""Per-variant: is the decline GLOBAL (purging) or CLIMATE-GRADED (hitchhiking)?

Restricted insertions decline in every garden while equally-restricted indels do not. That is a
population-level statement. This splits it per variant and asks which of two signatures each one
carries, using the residuals from `_sv_temporal_climate_matched.py` (s minus the median s of SNPs
matched on BOTH p0 and carrier-origin climate):

  GLOBAL          negative in significantly more of the 31 gardens than chance (sign test).
                  A deleterious insertion should lose ground everywhere, cold gardens included.
  CLIMATE-GRADED  its per-garden residual correlates with garden bio1 (Spearman), i.e. it falls
                  where warm-origin lineages win and rises where cold-origin ones do -- the
                  hitchhiking signature.

The two are not exclusive: a variant can be both (deleterious AND lineage-linked).

CALIBRATION. Every count is paired with the same count in the INDEL subsample, split the same way.
Indels are the reference class: restricted indels show no net decline, so their rate at any
threshold is the false-positive rate for these tests.

Env: kmate. Writes results/sv_adaptive/sv_purging_classify_{variants,summary}.csv.
"""
import os, sys
import numpy as np, pandas as pd
import scipy.sparse as sp
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
STORE = lib.AF_STORE
MIN_P0, SV_BP, NSUB_IND, ALPHA = 0.02, 50, 25000, 0.05
MIN_MAC = int(os.environ.get("SV_MIN_MAC", 12))
TAG = "" if MIN_MAC == 12 else f"_mac{MIN_MAC}"
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
UNIT = os.environ.get("SV_UNIT", "allele")         # "event": read the units scored by _sv_temporal_climate_matched.py
if UNIT == "event":
    FILTER = "record"; TAG = f"_event{MIN_MAC}"
if CRTAG:
    assert FILTER == "record", "SV_CALL_MIN only applies to the record / event filter"
    TAG += CRTAG


def record_state():
    """per-variant arrays for the record (allele) unit: af_store rows passing the filter."""
    F = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    grp = F.group.to_numpy()
    idx = np.load(f"{STORE}/index_nonsnp.npz")
    ch = idx["chrom"].astype("U5"); pos = idx["pos"].astype(np.int64)
    dl = idx["alt_len"].astype(np.int64) - idx["ref_len"].astype(np.int64)
    p0 = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
    keep = (panel_keep("nonsnp", ch, pos, MIN_MAC) & (p0 >= MIN_P0) & (p0 <= 1 - MIN_P0))
    cols = np.where(keep)[0]
    dln = dl[cols]; is_ins = dln > SV_BP; is_del = -dln > SV_BP
    is_ind = (np.abs(dln) >= 1) & (np.abs(dln) <= SV_BP)

    n = {g: np.zeros(cols.size, np.int64) for g in ("cold", "temperate", "warm")}
    tot = np.zeros(cols.size, np.int64)
    ncal = np.zeros(cols.size, np.int64)                      # founders genotyped (call rate)
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
            vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
            ncal[sel] = np.asarray((vc[:, nonsnp[cols[sel] - off]] > 0).sum(0)).ravel().astype(np.int64)
        off += nrow
        print(f"  [{C}] done", flush=True)

    ok = tot > 0
    rng = np.random.default_rng(0)
    ii = np.where(is_ind & ok)[0]
    i_sub = np.sort(rng.choice(ii, min(NSUB_IND, ii.size), replace=False))
    CH, POS, P0, AFROW = ch[cols], pos[cols], p0[cols], cols
    return_state = (cols, dln, is_ins, is_del, n, tot, ncal, i_sub, CH, POS, P0, AFROW)
    return return_state


def event_state():
    """the same per-variant arrays, for SV_UNIT=event: units as scored by the temporal step."""
    F = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    EU = np.load(f"{OUT}/nonsnp_event_units.npz")
    PS = np.load(f"{OUT}/sv_temporal_climate_matched_persite{TAG}.npz")
    cols = PS["cols_non"]; i_sub = PS["i_sub"]
    idx = np.load(f"{STORE}/index_nonsnp.npz")
    p0r = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
    A = sp.csr_matrix((np.ones(EU["indices"].size), EU["indices"], EU["indptr"]), shape=(EU["rep"].size, p0r.size))
    p0u = np.minimum(A[cols] @ p0r, 1.0)
    dln = EU["dlen"][cols]; is_ins = dln > SV_BP; is_del = -dln > SV_BP
    n = {"cold": EU["nc"][cols].astype(np.int64), "warm": EU["nw"][cols].astype(np.int64)}
    tot = EU["n"][cols].astype(np.int64); ncal = EU["ncal"][cols].astype(np.int64)
    assert np.array_equal(PS["p0_ins"], p0u[is_ins].astype(np.float32)), "unit order differs from temporal step"
    return (cols, dln, is_ins, is_del, n, tot, ncal, i_sub,
            EU["chrom"][cols].astype(str), EU["pos"][cols].astype(np.int64), p0u, EU["rep"][cols].astype(np.int64))


def main():
    F = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    (cols, dln, is_ins, is_del, n, tot, ncal, i_sub, CH, POS, P0, AFROW) = (
        event_state() if UNIT == "event" else record_state())

    def origin(midx):
        c, w = n["cold"][midx], n["warm"][midx]
        lab = np.full(midx.size, "other", dtype=object)
        lab[(c > 0) & (w > 0)] = "shared"
        lab[(c >= 2) & (w == 0)] = "cold-specific"
        lab[(w >= 2) & (c == 0)] = "warm-specific"
        return lab

    P = np.load(f"{OUT}/sv_temporal_climate_matched_persite{TAG}.npz")
    D = pd.read_csv(f"{OUT}/sv_temporal_climate_matched_sites{TAG}.csv").sort_values("bio1")
    sites = D.site.astype(int).tolist(); b1 = D.bio1.to_numpy()

    rows = []
    for cls, key, midx in (("insertion", "ins", np.where(is_ins)[0]),
                           ("deletion", "del", np.where(is_del)[0]),
                           ("indel", "indel", i_sub)):
        A = np.vstack([P[f"{s}_res_{key}_matched"] for s in sites])       # [gardens x variants]
        lab = origin(midx)
        okv = ~np.isnan(A)
        ngard = okv.sum(0)
        nneg = np.nansum(A < 0, axis=0)
        for j in range(A.shape[1]):
            if ngard[j] < 25:
                continue
            g = okv[:, j]
            # GLOBAL: negative in more gardens than chance
            pg = stats.binomtest(int(nneg[j]), int(ngard[j]), 0.5, alternative="greater").pvalue
            # CLIMATE-GRADED: residual tracks garden temperature
            r = stats.spearmanr(b1[g], A[g, j])
            lr = stats.linregress(b1[g], A[g, j])     # residual per degree C of garden bio1
            k_ = midx[j]
            rows.append(dict(cls=cls, origin=lab[j], col=int(k_), ridx=int(j), afrow=int(AFROW[k_]),
                             unit=int(cols[k_]) if UNIT == "event" else -1,
                             chrom=str(CH[k_]), pos=int(POS[k_]),
                             size_bp=int(abs(dln[k_])), n_carriers=int(tot[k_]),
                             n_cold=int(n["cold"][k_]), n_warm=int(n["warm"][k_]),
                             call_rate=ncal[k_] / len(F), p0=float(P0[k_]),
                             n_gardens=int(ngard[j]), frac_neg=nneg[j] / ngard[j],
                             mean_res=float(np.nanmean(A[:, j])), median_res=float(np.nanmedian(A[:, j])),
                             p_global=pg, rho_bio1=r.statistic, p_rho=r.pvalue,
                             slope_bio1=lr.slope, p_slope=lr.pvalue))
        print(f"  [{cls}] {A.shape[1]:,} variants scored", flush=True)

    V = pd.DataFrame(rows)
    V["is_global"] = V.p_global < ALPHA
    V["is_graded"] = (V.p_rho < ALPHA) & (V.rho_bio1 < 0)      # more purged at hot gardens
    V["call"] = np.where(V.is_global & V.is_graded, "both",
                 np.where(V.is_global, "global only",
                  np.where(V.is_graded, "climate-graded only", "neither")))
    V.to_csv(f"{OUT}/sv_purging_classify_variants{TAG}.csv", index=False)

    S = (V.groupby(["cls", "origin", "call"]).size().rename("n").reset_index())
    S.to_csv(f"{OUT}/sv_purging_classify_summary{TAG}.csv", index=False)
    print(f"\n{'class':<11}{'origin':<15}{'n':>6}{'global only':>13}{'graded only':>13}{'both':>8}{'neither':>9}")
    for cls in ("insertion", "deletion", "indel"):
        for o in ("cold-specific", "warm-specific", "shared"):
            d = V[(V.cls == cls) & (V.origin == o)]
            if len(d) < 20:
                continue
            c = d.call.value_counts()
            f = lambda k: f"{c.get(k,0)} ({100*c.get(k,0)/len(d):.0f}%)"
            print(f"{cls:<11}{o:<15}{len(d):>6}{f('global only'):>13}{f('climate-graded only'):>13}"
                  f"{f('both'):>8}{f('neither'):>9}")
        print()
    print(f"[wrote] {OUT}/sv_purging_classify_{{variants,summary}}{TAG}.csv")


if __name__ == "__main__":
    main()
