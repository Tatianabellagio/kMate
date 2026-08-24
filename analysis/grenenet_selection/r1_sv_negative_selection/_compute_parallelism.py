#!/usr/bin/env python
"""PARALLELISM across replicate plots (AF-vapeR / PicMin philosophy), per variant, SNP vs indel vs SV
(user 2026-07-03). Uses the ~10-12 replicate PLOTS within each site as parallel populations (same
founding pool + climate; each plot an independent pool-seq -> independent founder-h -> independent
allele-frequency-change). Consistency across plots = drift-controlled selection.

Per variant per site (K plots): per-plot logit-slope vector S[K].
  parallelism  rho = mean(S)^2 / mean(S^2)   in [0,1]  (0 = idiosyncratic/drift; 1 = fully parallel;
                                                        the rank-1 AF-vapeR eigenvalue analog)
  signed dir   = sign(mean(S));   among-plot z = mean/(sd/sqrt(K))
RESPONDER at a site = rho in the top RESP_Q within its p0-decile POOLED across classes (class-agnostic,
frequency-controlled). REPEATABILITY (PicMin-style) = # sites where the variant is a responder.

Saves per variant (all SV+indel-subsample+SNP-subsample): p0, class, isdel, n_sites, resp_count,
mean_rho, mean_absz, mean_signed_slope, and the per-site signed mean-slope matrix for a climate cross.
Env: kmate.  Writes analysis/grenenet_selection/sv_adaptive/parallelism.npz .
"""
import os, sys, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

os.chdir("/global/scratch/users/tbellg/kmate")
STORE = lib.AF_STORE; PM = f"{lib.GEA}/pool_matrices"
MIN_P0, SV_BP, MIN_MAC, NBIN, RESP_Q, EPS = 0.02, 50, 12, 10, 0.10, 1e-3
logit = lambda p: np.log(np.clip(p, EPS, 1-EPS)/(1-np.clip(p, EPS, 1-EPS)))


def read_rows(path, rows, cols):
    """Read specific ROWS of a [nrow x ncol] float32 .npy via seek+fromfile (sequential per row —
    ~150 MB/s vs ~1 MB/s for memmap page-faulted random access on this filesystem). -> [len(rows) x len(cols)]."""
    with open(path, "rb") as fh:
        v = np.lib.format.read_magic(fh)
        shp, fo, dt = np.lib.format.read_array_header_1_0(fh); off = fh.tell(); ncol = shp[1]
        out = np.empty((len(rows), cols.size), np.float64)
        for i, r in enumerate(rows):
            fh.seek(off + int(r) * ncol * dt.itemsize)
            out[i] = np.fromfile(fh, dtype=dt, count=ncol)[cols]
    return out


def plot_slopes(site, kind, cols, p0c):
    """S [K x ncols] per-plot logit-slopes across the site's K plots (gen0 = shared p0c)."""
    L0 = logit(p0c); plots = {}
    for g in (1, 2, 3):
        try: meta = pd.read_csv(f"{PM}/pool_gen{g}_{kind}.meta.csv")
        except FileNotFoundError: continue
        m4 = meta[meta.site == site]
        if len(m4) == 0: continue
        M = read_rows(f"{PM}/pool_gen{g}_{kind}_af.npy", m4.index.to_numpy(), cols)
        for j, plot in enumerate(m4["plot"].to_numpy()):
            plots.setdefault(int(plot), []).append((float(g), logit(M[j])))
    S = []
    for plot, pts in plots.items():
        ts = np.array([0.0] + [t for t, _ in pts])
        if ts.size < 2: continue
        Y = np.vstack([L0] + [y for _, y in pts]); tc = ts - ts.mean()
        S.append((tc @ Y) / (tc ** 2).sum())
    return np.vstack(S)


def main():
    idx_snp = np.load(f"{STORE}/index_snp.npz"); idx_non = np.load(f"{STORE}/index_nonsnp.npz")
    ch_non = idx_non["chrom"].astype("U5"); pos_non = idx_non["pos"].astype(np.int64)
    dlen = np.abs(idx_non["alt_len"].astype(np.int64) - idx_non["ref_len"].astype(np.int64))
    is_del = idx_non["ref_len"].astype(np.int64) > idx_non["alt_len"].astype(np.int64)
    ch_snp = idx_snp["chrom"].astype("U5"); pos_snp = idx_snp["pos"].astype(np.int64)
    common_non = lib.founder_panel_keep(ch_non, pos_non, min_mac=MIN_MAC)
    common_snp = lib.founder_panel_keep(ch_snp, pos_snp, min_mac=MIN_MAC)
    p0_non = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
    p0_snp = np.load(f"{STORE}/p0_snp.npy").astype(np.float64)
    snp_m = common_snp & (p0_snp >= MIN_P0) & (p0_snp <= 1 - MIN_P0)
    non_m = (dlen >= 1) & common_non & (p0_non >= MIN_P0) & (p0_non <= 1 - MIN_P0)
    rng = np.random.default_rng(0)
    # subsample SNP + indel columns (keep all SV) to bound the per-site output matrices
    cols_sv = np.where(non_m & (dlen > SV_BP))[0]
    ci = np.where(non_m & (dlen >= 1) & (dlen <= SV_BP))[0]; cols_in = rng.choice(ci, 50000, replace=False)
    cs = np.where(snp_m)[0]; cols_sn = rng.choice(cs, 50000, replace=False)
    # unified nonsnp column set = SV + indel-subsample (one pass over nonsnp matrices)
    cols_non = np.concatenate([cols_sv, cols_in]); is_sv = np.concatenate([np.ones(cols_sv.size, bool), np.zeros(cols_in.size, bool)])
    isdel_non = is_del[cols_non]; p0_non_k = p0_non[cols_non]; p0_sn_k = p0_snp[cols_sn]
    print(f"cols: SV={cols_sv.size:,} indel(sub)={cols_in.size:,} SNP(sub)={cols_sn.size:,}", flush=True)

    clim = pd.concat([pd.read_csv(f"{PM}/pool_gen{g}_snp.meta.csv") for g in (1,2,3)],
                     ignore_index=True).groupby("site")[["bio1", "bio18"]].mean()
    sites = sorted({int(s) for f in glob.glob(f"{PM}/pool_gen*_snp.meta.csv") for s in pd.read_csv(f).site.unique()})
    # global p0 deciles from SNP for the responder threshold binning
    p0q = np.quantile(p0_snp[snp_m], np.linspace(0, 1, NBIN + 1)); p0q[-1] += 1e-6
    binf = lambda p: np.clip(np.digitize(p, p0q[1:-1]), 0, NBIN - 1)
    b_non = binf(p0_non_k); b_sn = binf(p0_sn_k)

    n_non = cols_non.size; n_sn = cols_sn.size
    rho_sum = {"non": np.zeros(n_non), "sn": np.zeros(n_sn)}
    absz_sum = {"non": np.zeros(n_non), "sn": np.zeros(n_sn)}
    resp = {"non": np.zeros(n_non), "sn": np.zeros(n_sn)}
    nsite = {"non": np.zeros(n_non), "sn": np.zeros(n_sn)}
    signed = {"non": [], "sn": []}; perz = {"non": [], "sn": []}; perrho = {"non": [], "sn": []}
    site_bio1 = []; site_bio18 = []; used_sites = []

    def stats(S):
        K = S.shape[0]; m = S.mean(0); sd = S.std(0, ddof=1)
        rho = m ** 2 / np.maximum((S ** 2).mean(0), 1e-12)
        z = np.divide(m, sd / np.sqrt(K), out=np.zeros_like(m), where=sd > 0)
        return m, rho, np.abs(z)

    for site in sites:
        if site not in clim.index: continue
        try:
            S_sn = plot_slopes(site, "snp", cols_sn, p0_sn_k)
            S_no = plot_slopes(site, "nonsnp", cols_non, p0_non_k)
        except Exception as e:
            print(f"  site {site}: skip ({e})"); continue
        if S_sn.shape[0] < 3: continue                    # need >=3 plots for parallelism
        m_sn, rho_sn, az_sn = stats(S_sn); m_no, rho_no, az_no = stats(S_no)
        # responder threshold per p0-bin from POOLED rho (class-agnostic)
        for b in range(NBIN):
            pool = np.concatenate([rho_sn[b_sn == b], rho_no[b_non == b]])
            if pool.size < 20: continue
            thr = np.quantile(pool, 1 - RESP_Q)
            resp["sn"][b_sn == b] += (rho_sn[b_sn == b] > thr)
            resp["non"][b_non == b] += (rho_no[b_non == b] > thr)
        for tag, rho, az, m, nk in (("sn", rho_sn, az_sn, m_sn, n_sn), ("non", rho_no, az_no, m_no, n_non)):
            rho_sum[tag] += rho; absz_sum[tag] += az; nsite[tag] += 1; signed[tag].append(m.astype(np.float32))
            perz[tag].append(az.astype(np.float32)); perrho[tag].append(rho.astype(np.float32))
        site_bio1.append(float(clim.loc[site, "bio1"])); site_bio18.append(float(clim.loc[site, "bio18"])); used_sites.append(site)
        print(f"  site {site:>2} plots={S_sn.shape[0]} done", flush=True)

    out = dict(
        p0_sv=p0_non_k[is_sv].astype(np.float32), isdel_sv=isdel_non[is_sv],
        resp_sv=resp["non"][is_sv], nsite_sv=nsite["non"][is_sv], rho_sv=(rho_sum["non"]/nsite["non"])[is_sv].astype(np.float32),
        absz_sv=(absz_sum["non"]/nsite["non"])[is_sv].astype(np.float32),
        p0_indel=p0_non_k[~is_sv].astype(np.float32),
        resp_indel=resp["non"][~is_sv], nsite_indel=nsite["non"][~is_sv], rho_indel=(rho_sum["non"]/nsite["non"])[~is_sv].astype(np.float32),
        p0_snp=p0_sn_k.astype(np.float32), resp_snp=resp["sn"], nsite_snp=nsite["sn"],
        rho_snp=(rho_sum["sn"]/nsite["sn"]).astype(np.float32), absz_snp=(absz_sum["sn"]/nsite["sn"]).astype(np.float32),
        site_bio1=np.array(site_bio1), site_bio18=np.array(site_bio18), used_sites=np.array(used_sites),
        signed_sv=np.vstack(signed["non"]).T[is_sv].astype(np.float32),      # [n_sv x nsite] for climate cross
        signed_snp=np.vstack(signed["sn"]).T.astype(np.float32),
        # per-site drift-controlled parallelism (|z| and rho) per variant -> PicMin input [nvar x nsite]
        z_sv=np.vstack(perz["non"]).T[is_sv].astype(np.float32), z_indel=np.vstack(perz["non"]).T[~is_sv].astype(np.float32),
        z_snp=np.vstack(perz["sn"]).T.astype(np.float32),
        rhops_sv=np.vstack(perrho["non"]).T[is_sv].astype(np.float32), rhops_snp=np.vstack(perrho["sn"]).T.astype(np.float32))
    np.savez_compressed(f"{lib.GEA}/sv_adaptive/parallelism.npz", **out)

    from scipy import stats as st
    def matchcmp(a_val, a_p0, b_val, b_p0):
        e = np.quantile(np.concatenate([a_p0, b_p0]), np.linspace(0, 1, 26)); e[-1] += 1e-9
        ab = np.digitize(a_p0, e[1:-1]); bb = np.digitize(b_p0, e[1:-1]); m = []
        for k in range(25):
            pool = b_val[bb == k]; nn = int(round((ab == k).mean() * 20000))
            if nn and pool.size: m.append(rng.choice(pool, nn, replace=True))
        return np.concatenate(m)
    msnp_resp = matchcmp(out["resp_sv"], out["p0_sv"], out["resp_snp"], out["p0_snp"])
    print("\n=== PARALLEL-RESPONDER RATE: SV vs frequency-matched SNP ===")
    print(f"  mean responder-count / #sites: SV={out['resp_sv'].mean()/np.median(out['nsite_sv']):.3f}  "
          f"matched-SNP={msnp_resp.mean()/np.median(out['nsite_snp']):.3f}")
    print(f"  mean rho: SV={np.mean(out['rho_sv']):.3f} indel={np.mean(out['rho_indel']):.3f} SNP={np.mean(out['rho_snp']):.3f}")
    print(f"  frac 'repeatable' (responder in >=1/3 of sites): SV={np.mean(out['resp_sv']>=out['nsite_sv']/3):.3f} "
          f"matched-SNP={np.mean(msnp_resp>=np.median(out['nsite_snp'])/3):.3f}")
    print("[wrote] parallelism.npz")


if __name__ == "__main__":
    main()
