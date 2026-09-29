#!/usr/bin/env python
"""Per-garden allele frequencies for the GrENE-Net panels (Fig. 2c-d, Fig. 3), straight from kMate.

ONE data source, no recomputation: kMate's per-pool allele frequencies as stored in
`common/results/pool_matrices/pool_gen9_{nonsnp,snp}_af.npy` -- "gen 9" = each plot's LAST sampled
generation (352 pools, all 31 gardens) -- and the founding frequency p0 = mean of the 8 seed-mix kMate
runs (`af_store/p0_{nonsnp,snp}.npy`). A garden's frequency of a variant = mean over its plots of the
pool frequency (plots equal-weighted; NaN pools dropped).

Units: SV EVENT = Truvari-collapsed alleles (`nonsnp_event_units.npz`, same order as
`sv_event_catalog.csv`); event frequency in a pool = sum of its member alleles' frequencies, capped at
1 (NaN if any member is NaN). Short indels per record (random 10%); SNPs per record (random 10%).
Filters as everywhere: called in >= 90% of the 231 founders, >= 2 carriers. Labels as in the catalog:
context (first match wins), origin (shared / cold-specific / warm-specific / other), carrier count n,
mean carrier-origin bio1 (b1).

GARDEN_STAT=s: the same table with, in cols 1..31, the per-generation selection coefficient instead of the
frequency -- per plot s = [logit(p_T) - logit(p_0)] / T (T = the plot's last sampled generation, clip 1e-3),
averaged over a garden's plots; col 0 stays p0. Written to garden_s_event.npz (garden_af_event.npz untouched).
GARDEN_STAT=slope: s as the SLOPE ACROSS GENERATIONS -- per plot, OLS slope of logit allele frequency on
generation (t = 0 seed mix + every sampled generation 1-3, `pool_gen{1,2,3}_*`), averaged over a garden's
plots. Written to garden_sslope_event.npz.
GARDEN_STAT=gslope (the one to use): the same slope at the GARDEN level -- per generation the garden's frequency is the
mean over its plots, THEN logit and OLS slope on generation. Plot-level frequencies are noisy and often exactly 0, and
the logit of those dominates a plot-then-average estimate. Written to garden_sgslope_event.npz.

Env: kmate. Writes results/sv_adaptive/garden_af_event.npz
  p [variants x 32]  col 0 = p0, cols 1..31 = gardens (sites sorted by bio1 as in `sites`)
  kind, ctx, n, nc, nw, b1, origin, seq_class, cat_row (catalog row for events, -1 otherwise)
  sites, bio1, zone
and sv_event_catalog_p_poolaf.npy: the same [catalog rows x 32] for every passing event (NaN otherwise).
"""
import os, sys, importlib.util
import numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
_H = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location("sfsctx", os.path.join(_H, "_sv_sfs_context.py"))
sc = importlib.util.module_from_spec(_s); _s.loader.exec_module(sc)

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
ST = lib.AF_STORE; PM = f"{lib.GEA}/common/results/pool_matrices"
N, CALL_MIN, MIN_CAR, SV_BP, KEEP, SEED = 231, 0.9, 2, 50, 0.10, 0
CTXI = {c: i for i, c in enumerate(sc.CTX)}
STAT = os.environ.get("GARDEN_STAT", "af")          # "af" (default) | "s": endpoint s | "slope": s as the slope across generations
NO_T0 = STAT == "gslope_no0"                        # gslope without the seed-mix point: p0 enters nowhere; needs >= 2 sampled generations
if NO_T0: STAT = "gslope"
EPS_S = 1e-3                                        # frequency clip before the logit
UNIT = os.environ.get("SV_UNIT", "event")           # "event" (Truvari-collapsed, default) | "allele": every SV record on its own


def pick_cols(path, cols, block=8):
    """X[:, cols] from a big [pools x records] .npy without loading it: sequential row blocks read with
    plain file reads (mmap page faults crawl on Lustre; read() streams at disk speed). Only the picked
    columns are kept in RAM (this allocation has 16 GB; the SNP matrix alone is 8.8 GB)."""
    with open(path, "rb") as f:
        version = np.lib.format.read_magic(f)
        reader = {(1, 0): np.lib.format.read_array_header_1_0, (2, 0): np.lib.format.read_array_header_2_0}[version]
        shape, fortran, dtype = reader(f); assert not fortran
        nrow, ncol = shape; out = np.empty((nrow, cols.size), np.float32); off0 = f.tell()
        for i in range(0, nrow, block):
            nb = min(block, nrow - i); f.seek(off0 + i * ncol * dtype.itemsize)
            buf = np.fromfile(f, dtype=dtype, count=nb * ncol).reshape(nb, ncol)
            out[i:i + nb] = buf[:, cols]
    return out


def site_means(X, site_of_pool, sites):
    """mean over a site's pools, NaN-aware. X: [pools x variants] -> [variants x n_sites]"""
    out = np.full((X.shape[1], len(sites)), np.nan, np.float32)
    for j, s in enumerate(sites):
        m = site_of_pool == s
        if m.any():
            out[:, j] = np.nanmean(X[m], 0)
    return out


def to_stat(X, p0, gen):
    """STAT=af: the pool frequencies themselves. STAT=s: per pool, the haploid (selfing) selection
    coefficient per generation,  s = [logit(p_T) - logit(p_0)] / T,  T = that plot's last sampled
    generation (1-3), frequencies clipped to [EPS_S, 1 - EPS_S]; gardens then average s over plots."""
    if STAT != "s":
        return X
    lg = lambda p: np.log(np.clip(p, EPS_S, 1 - EPS_S) / (1 - np.clip(p, EPS_S, 1 - EPS_S)))
    return (lg(np.asarray(X, np.float64)) - lg(np.asarray(p0, np.float64))[None, :]) / gen[:, None]


def lgt(p):
    p = np.clip(np.asarray(p, np.float64), EPS_S, 1 - EPS_S); return np.log(p / (1 - p))


def slopes(pool_by_gen, metas, p0):
    """STAT=slope. Per plot and variant: OLS slope of logit allele frequency on generation, t = 0 (seed mix, p0) plus
    every generation the plot was sampled in (1-3). pool_by_gen {g: [pools_g x variants]}, metas {g: meta with site, plot}.
    Returns [plots x variants] slopes and the site of each plot (plots with >= 1 sampled generation)."""
    L0 = lgt(p0); plots = {}
    for g, meta in metas.items():
        for r, (s, pl) in enumerate(zip(meta["site"].to_numpy(), meta["plot"].to_numpy())):
            plots.setdefault((s, pl), []).append((g, r))
    out = np.empty((len(plots), L0.size), np.float32); site_of = np.empty(len(plots), np.int64)
    for i, ((s, pl), pts) in enumerate(plots.items()):
        ts = np.array([0.0] + [g for g, _ in pts]); tc = ts - ts.mean()
        Y = np.vstack([L0] + [lgt(pool_by_gen[g][r]) for g, r in pts])
        out[i] = (tc @ Y) / (tc ** 2).sum(); site_of[i] = s
    return out, site_of


def garden_slopes(pool_by_gen, metas, p0, sites, use_t0=True):
    """STAT=gslope. s as the slope across generations AT THE GARDEN LEVEL: per generation, a garden's frequency = mean
    over its plots sampled in that generation (plots are replicates of one population; averaging BEFORE the logit keeps
    plot-level zeros and sampling noise from dominating); then, per variant, the OLS slope of logit frequency on
    generation through t = 0 (seed mix) and every generation the garden was sampled in. Returns [variants x n_sites]."""
    L0 = lgt(p0); out = np.full((L0.size, len(sites)), np.nan, np.float32)
    for j, s in enumerate(sites):
        ts, Y = ([0.0], [L0]) if use_t0 else ([], [])
        for g, meta in metas.items():
            m = meta["site"].to_numpy() == s
            if m.any():
                ts.append(float(g)); Y.append(lgt(np.nanmean(pool_by_gen[g][m], 0)))
        if len(ts) > 1:
            tc = np.array(ts) - np.mean(ts); out[:, j] = (tc @ np.vstack(Y)) / (tc ** 2).sum()
    return out


def main():
    rng = np.random.default_rng(SEED)
    FO = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    GC = (FO.group == "cold").to_numpy(); GW = (FO.group == "warm").to_numpy(); B1 = FO.bio1.to_numpy()
    meta_n = pd.read_csv(f"{PM}/pool_gen9_nonsnp.meta.csv"); meta_s = pd.read_csv(f"{PM}/pool_gen9_snp.meta.csv")
    assert meta_n.pool.tolist() == meta_s.pool.tolist(), "non-SNP and SNP pool matrices differ in pool order"
    site_bio1 = meta_n.groupby("site").bio1.first().sort_values(); sites = site_bio1.index.to_numpy(); bio1 = site_bio1.to_numpy()
    zone = pd.qcut(bio1, 3, labels=["cold", "mid", "hot"]).astype(str).to_numpy()
    sop = meta_n.site.to_numpy(); gen = meta_n.generation.to_numpy(np.float64)
    print(f"gen-9 pools: {len(meta_n)} over {len(sites)} gardens (each plot at its last sampled generation)", flush=True)
    F_n = f"{PM}/pool_gen9_nonsnp_af.npy"; F_s = f"{PM}/pool_gen9_snp_af.npy"
    p0n = np.load(f"{ST}/p0_nonsnp.npy").astype(np.float64); p0s = np.load(f"{ST}/p0_snp.npy").astype(np.float64)
    F = {k: [] for k in ("kind", "ctx", "n", "nc", "nw", "b1", "origin", "seq_class", "cat_row", "p")}

    def put(kind, ctx, n, nc, nw, b1, origin, seq_class, cat_row, p):
        F["kind"].append(np.full(n.size, kind)); F["ctx"].append(ctx.astype(np.int8)); F["n"].append(n.astype(np.int16))
        F["nc"].append(nc.astype(np.int16)); F["nw"].append(nw.astype(np.int16)); F["b1"].append(b1.astype(np.float32))
        F["origin"].append(np.asarray(origin).astype("U24")); F["seq_class"].append(np.asarray(seq_class).astype("U24"))   # fixed-width: no pickle, readable by any numpy
        F["cat_row"].append(cat_row.astype(np.int64)); F["p"].append(p.astype(np.float32))

    # ---- SV events (catalog order == event-unit order; verified on n_carriers / n_called) ----
    D = pd.read_csv(f"{OUT}/sv_event_catalog.csv"); EU = np.load(f"{OUT}/nonsnp_event_units.npz")
    dlen = EU["dlen"]; svu = np.where(np.abs(dlen) > SV_BP)[0]
    assert svu.size == len(D) and np.array_equal(EU["n"][svu], D.n_carriers.to_numpy()) and np.array_equal(EU["ncal"][svu], D.n_called.to_numpy())
    rows_pass = np.where(D.passes.to_numpy())[0]; units = svu[rows_pass]
    indptr, indices = EU["indptr"], EU["indices"]
    if UNIT == "allele":
        # every SV record (|dlen| > 50) as its own variant: own carriers >= 2, own call rate >= 90%, own context / origin
        idx_n = np.load(f"{ST}/index_nonsnp.npz"); ch_n = idx_n["chrom"].astype("U5"); IVa = sc.layers(); AL = []; offn = 0
        for C in lib.CHROMS:
            base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
            meta = np.load(f"{base}.meta.npz", allow_pickle=True)
            pos = meta["pos"].astype(np.int64); rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
            nons = np.where(~((rl == 1) & (al == 1)))[0]; nC = int((ch_n == C).sum()); assert nons.size == nC
            sel = np.where(np.abs(al[nons] - rl[nons]) > SV_BP)[0]; cols = nons[sel]
            vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
            P = (vp[:, cols] > 0).tocsc(); c = np.asarray(P.sum(0)).ravel(); k = np.asarray((vc[:, cols] > 0).sum(0)).ravel()
            ok = np.where((k >= CALL_MIN * N) & (c >= MIN_CAR) & (c <= N - MIN_CAR))[0]; P = P[:, ok]; cols = cols[ok]
            nc_ = np.asarray(P[GC].sum(0)).ravel(); nw_ = np.asarray(P[GW].sum(0)).ravel()
            AL.append(dict(row=offn + sel[ok], kind=np.where(al[cols] > rl[cols], "ins", "del"), n=c[ok], nc=nc_, nw=nw_,
                           b1=np.asarray(P.astype(np.float64).T @ B1).ravel() / c[ok],
                           ctx=sc.context(IVa, np.full(cols.size, C, dtype=object), pos[cols], pos[cols] + np.maximum(rl[cols] - 1, 0))))
            offn += nC
        AL = {k_: np.concatenate([a[k_] for a in AL]) for k_ in AL[0]}
        mem = [np.array([r]) for r in AL["row"]]; need = np.unique(AL["row"])
        print(f"[alleles] {AL['row'].size:,} SV records pass on their own carriers / call rate", flush=True)
    else:
        mem = [indices[indptr[u]:indptr[u + 1]] for u in units]; need = np.unique(np.concatenate(mem))
    col = {r: j for j, r in enumerate(need)}
    M = sp.csr_matrix((np.ones(sum(len(m) for m in mem)), np.concatenate([[col[r] for r in m] for m in mem]),
                       np.r_[0, np.cumsum([len(m) for m in mem])]), shape=(len(mem), need.size))
    # ---- short indels: event units with 1 <= |dlen| <= 50 (single records), random 10% ----
    iu = np.where((np.abs(dlen) >= 1) & (np.abs(dlen) <= SV_BP) & (EU["n"] >= MIN_CAR) & (EU["n"] <= N - MIN_CAR) & (EU["ncal"] >= CALL_MIN * N))[0]
    iu = np.sort(rng.choice(iu, int(KEEP * iu.size), replace=False))
    assert np.all(indptr[iu + 1] - indptr[iu] == 1); irow = indices[indptr[iu]]
    allrows = np.unique(np.r_[need, irow]); pos_in = {r: j for j, r in enumerate(allrows)}
    XN = pick_cols(F_n, allrows); print(f"  non-SNP matrix: one pass, {allrows.size:,} columns kept", flush=True)
    def events_of(XN_):
        X = XN_[:, [pos_in[r] for r in need]].astype(np.float64)                 # pools x member rows
        bad = np.isnan(X); Xz = np.where(bad, 0.0, X)
        ev = np.minimum(np.asarray((M @ Xz.T)), 1.0).T                            # pools x events, capped sum
        ev[np.asarray((M @ bad.T.astype(np.float64))).T > 0] = np.nan             # any member NaN -> NaN
        return ev
    ev_pool = events_of(XN)
    ev_p0 = np.minimum(M @ p0n[need], 1.0)
    GENS = (1, 2, 3); XNG = {}
    if STAT in ("slope", "gslope"):                                               # per-generation pool frequencies
        metas = {g: pd.read_csv(f"{PM}/pool_gen{g}_nonsnp.meta.csv") for g in GENS}
        for g in GENS:
            assert metas[g].pool.tolist() == pd.read_csv(f"{PM}/pool_gen{g}_snp.meta.csv").pool.tolist()
            XNG[g] = pick_cols(f"{PM}/pool_gen{g}_nonsnp_af.npy", allrows); print(f"  gen {g}: {len(metas[g])} pools", flush=True)
        if STAT == "gslope":
            P_ev = np.c_[ev_p0, garden_slopes({g: events_of(XNG[g]) for g in GENS}, metas, ev_p0, sites, use_t0=not NO_T0)]
        else:
            S_, sop_pl = slopes({g: events_of(XNG[g]) for g in GENS}, metas, ev_p0)
            P_ev = np.c_[ev_p0, site_means(S_, sop_pl, sites)]
    else:
        P_ev = np.c_[ev_p0, site_means(to_stat(ev_pool, ev_p0, gen), sop, sites)]
    if UNIT == "allele":
        org = np.where((AL["nc"] > 0) & (AL["nw"] > 0), "shared", np.where((AL["nc"] >= 2) & (AL["nw"] == 0), "cold-specific",
                       np.where((AL["nw"] >= 2) & (AL["nc"] == 0), "warm-specific", "other")))
        put(AL["kind"], AL["ctx"], AL["n"], AL["nc"], AL["nw"], AL["b1"], org, np.full(AL["n"].size, ""), np.full(AL["n"].size, -1), P_ev)
    else:
        d = D.iloc[rows_pass]
        put(d.kind.to_numpy(), d.context.map(CTXI).to_numpy(), d.n_carriers.to_numpy(), d.n_cold.to_numpy(), d.n_warm.to_numpy(),
            d.carrier_bio1.to_numpy(), d.origin.to_numpy(), d.seq_class.fillna("").to_numpy(), rows_pass, P_ev)
        PC = np.full((len(D), 32), np.nan, np.float32); PC[rows_pass] = P_ev
    if STAT == "af" and UNIT == "event":
        np.save(f"{OUT}/sv_event_catalog_p_poolaf.npy", PC)
    print(f"[events] {rows_pass.size:,} passing events from {need.size:,} allele rows; p0 NaN {np.isnan(ev_p0).sum()}, garden NaN share {np.isnan(P_ev[:, 1:]).mean():.4f}", flush=True)

    IV = sc.layers()
    ictx = sc.context(IV, EU["chrom"][iu].astype(object), EU["pos"][iu].astype(np.int64), EU["pos"][iu].astype(np.int64))
    Xi = XN[:, [pos_in[r] for r in irow]]; del XN
    nc, nw, n = EU["nc"][iu], EU["nw"][iu], EU["n"][iu]
    origin = np.where((nc > 0) & (nw > 0), "shared", np.where((nc >= 2) & (nw == 0), "cold-specific", np.where((nw >= 2) & (nc == 0), "warm-specific", "other")))
    if STAT == "gslope":
        icol = [pos_in[r] for r in irow]; Pi = np.c_[p0n[irow], garden_slopes({g: XNG[g][:, icol] for g in GENS}, metas, p0n[irow], sites, use_t0=not NO_T0)]; del XNG
    elif STAT == "slope":
        icol = [pos_in[r] for r in irow]; S_, _ = slopes({g: XNG[g][:, icol] for g in GENS}, metas, p0n[irow]); Pi = np.c_[p0n[irow], site_means(S_, sop_pl, sites)]; del XNG
    else:
        Pi = np.c_[p0n[irow], site_means(to_stat(Xi, p0n[irow], gen), sop, sites)]
    put(np.full(iu.size, "indel"), ictx, n, nc, nw, EU["b1"][iu], origin, np.full(iu.size, ""), np.full(iu.size, -1), Pi)
    print(f"[indels] {iu.size:,} records (10% of passing)", flush=True)

    # ---- SNPs: per chromosome, af_store SNP rows == panel SNP columns in order (asserted in lib) ----
    idx_s = np.load(f"{ST}/index_snp.npz"); ch_s = idx_s["chrom"].astype("U5"); off = 0; SN = []
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64); rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        cols = np.where((rl == 1) & (al == 1))[0]; nC = int((ch_s == C).sum()); assert cols.size == nC
        sub = np.where(rng.random(cols.size) < KEEP)[0]                        # 10% of SNP records
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc(); vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        P = (vp[:, cols[sub]] > 0).tocsc(); Q = (vc[:, cols[sub]] > 0).tocsc()
        c = np.asarray(P.sum(0)).ravel(); k = np.asarray(Q.sum(0)).ravel()
        ok = np.where((k >= CALL_MIN * N) & (c >= MIN_CAR) & (c <= N - MIN_CAR))[0]
        P = P[:, ok]; sub = sub[ok]; c = c[ok]
        nc = np.asarray(P[GC].sum(0)).ravel(); nw = np.asarray(P[GW].sum(0)).ravel()
        b1 = np.asarray(P.astype(np.float64).T @ B1).ravel() / np.maximum(c, 1)
        ctx = sc.context(IV, np.full(sub.size, C, dtype=object), pos[cols[sub]], pos[cols[sub]])
        srow = off + sub
        origin = np.where((nc > 0) & (nw > 0), "shared", np.where((nc >= 2) & (nw == 0), "cold-specific", np.where((nw >= 2) & (nc == 0), "warm-specific", "other")))
        SN.append((srow, ctx, c, nc, nw, b1, origin)); print(f"[{C}] SNPs {sub.size:,} (10% of passing)", flush=True); off += nC
    srow = np.concatenate([t[0] for t in SN]); XS = pick_cols(F_s, srow) if STAT not in ("slope", "gslope") else None; print(f"  SNP matrix: {srow.size:,} columns kept", flush=True)
    cat = lambda i: np.concatenate([t[i] for t in SN])
    if STAT == "gslope":
        del XS; Ps = np.c_[p0s[srow], garden_slopes({g: pick_cols(f"{PM}/pool_gen{g}_snp_af.npy", srow) for g in GENS}, metas, p0s[srow], sites, use_t0=not NO_T0)]
    elif STAT == "slope":
        del XS; S_, _ = slopes({g: pick_cols(f"{PM}/pool_gen{g}_snp_af.npy", srow) for g in GENS}, metas, p0s[srow]); Ps = np.c_[p0s[srow], site_means(S_, sop_pl, sites)]
    else:
        Ps = np.c_[p0s[srow], site_means(to_stat(XS, p0s[srow], gen), sop, sites)]
    put(np.full(srow.size, "snp"), cat(1), cat(2), cat(3), cat(4), cat(5), cat(6), np.full(srow.size, ""), np.full(srow.size, -1), Ps)

    S = {k: np.concatenate(v) for k, v in F.items()}
    S["kind"] = S["kind"].astype("U8")
    fn = {"af": "garden_af_event.npz" if UNIT == "event" else "garden_af_allele.npz", "s": "garden_s_event.npz", "slope": "garden_sslope_event.npz", "gslope": "garden_sgslope_event.npz" if UNIT == "event" else "garden_sgslope_allele.npz"}[STAT]   # same rows, same order
    if NO_T0: fn = fn.replace("sgslope", "sgslope_no0")
    np.savez_compressed(f"{OUT}/{fn}", sites=sites, bio1=bio1, zone=zone.astype("U8"), **S)
    print("\nclass counts:"); print(pd.Series(S["kind"]).value_counts().to_string())
    print(f"[wrote] {OUT}/{fn}" + (" and sv_event_catalog_p_poolaf.npy" if STAT == "af" else ""))


if __name__ == "__main__":
    main()
