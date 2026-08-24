#!/usr/bin/env python
"""DECISIVE founder-level test of the "selection on SVs" conclusion (audit 2026-07-06).

In GLOBAL mode, per-variant pool AF is AF_v = h.V_pa[:,v]/h.V_called[:,v]: a deterministic
projection of the per-chrom founder mixture h through the FIXED founder x variant matrices.
So any per-variant statistic (the climate-slope beta from _compute_s_climate_slope.py) can
carry information ONLY through (i) which founders carry v and (ii) how those founders respond
to climate. It contains no per-SV selection signal by construction.

This script makes that concrete. Using ONLY:
  - gamma_f = d s_f / d climate  = per-FOUNDER climate response (slope of the founder selection
    coefficient S[site,f] on climate across the 31 sites; from selection_s_matrix.npz), and
  - the carrier sets from the founder panel var_pa (NO per-variant temporal data at all),
we build a founder-composition predictor for every panel variant:
      gbar_v = sum_{f in carriers} p0_f * gamma_f / sum_{f in carriers} p0_f
and ask:
  (T1) Across founders: do founders that carry MORE insertions have more negative gamma
       (decline in hot/arid)?  -> the founder-level insertion-load x climate pattern.
  (T2) Does gbar (founder composition ALONE) reproduce the SV-insertion-vs-SNP shift and the
       insertion/deletion asymmetry that the temporal beta shows?
  (T3) Does gbar predict the actual per-variant temporal beta (Spearman)? and does conditioning
       on gbar remove the ins/del difference in beta?  -> if yes, the "per-SV selection" is
       fully a founder-carrier-composition effect (audit finding M1/M4).

Env: kmate. Reads selection_s_matrix.npz, s_climate_slope.npz, panel/arch3, af_store.
Writes analysis/grenenet_gea/sv_adaptive/founder_load_test.npz + prints summary.
"""
import os, sys, glob
import numpy as np, pandas as pd
import scipy.sparse as sp
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

os.chdir("/global/scratch/users/tbellg/kmate")
STORE = lib.AF_STORE; GEA = lib.GEA; PROJ = "/global/scratch/users/tbellg/kmate"
PM = f"{GEA}/pool_matrices"
MIN_P0, SV_BP, MIN_MAC, NBIN, NDRAW = 0.02, 50, 12, 25, 20000
rng = np.random.default_rng(0)


def match_to(tp0, ss, sp0):
    """frequency-match ss (with freqs sp0) to the target freq dist tp0; return matched sample."""
    e = np.quantile(np.concatenate([tp0, sp0]), np.linspace(0, 1, NBIN + 1)); e[-1] += 1e-9
    tb = np.digitize(tp0, e[1:-1]); sb = np.digitize(sp0, e[1:-1]); out = []
    for b in range(NBIN):
        k = int(round((tb == b).mean() * NDRAW)); pool = ss[sb == b]
        if k > 0 and pool.size:
            out.append(rng.choice(pool, k, replace=True))
    return np.concatenate(out)


def main():
    # ---------- founder climate-response gamma_f (bio1, bio18) ----------
    Z = np.load(f"{GEA}/varexp/selection_s_matrix.npz", allow_pickle=True)
    S = Z["S"]                      # [31 sites x 231 founders]  founder selection coef per site
    founders = Z["founders"].astype("U6")
    sites = Z["sites"].astype(int)
    analyz = Z["analyzable"].astype(bool)
    p0f = Z["p0"].astype(np.float64)          # founder baseline freq (seedmix)
    # per-site climate (same source the temporal script used)
    cmeta = pd.concat([pd.read_csv(f"{PM}/pool_gen{g}_snp.meta.csv") for g in (1, 2, 3)],
                      ignore_index=True).groupby("site")[["bio1", "bio18"]].mean()
    clim = {k: np.array([cmeta.loc[s, k] if s in cmeta.index else np.nan for s in sites])
            for k in ("bio1", "bio18")}

    def founder_gamma(cvar):
        c = clim[cvar]; ok_site = np.isfinite(c)
        cc = c[ok_site] - c[ok_site].mean(); denom = float((cc ** 2).sum())
        g = np.full(founders.size, np.nan)
        for f in range(founders.size):
            col = S[ok_site, f]
            m = np.isfinite(col)
            if m.sum() < 5:
                continue
            cx = c[ok_site][m] - c[ok_site][m].mean()
            g[f] = float((cx @ col[m]) / (cx @ cx))
        return g
    gamma = {k: founder_gamma(k) for k in ("bio1", "bio18")}
    for k in ("bio1", "bio18"):
        gg = gamma[k][analyz]
        print(f"[gamma] {k}: {np.isfinite(gg).sum()} analyzable founders, "
              f"median={np.nanmedian(gg):+.4f} sd={np.nanstd(gg):.4f}")

    # ---------- per-chrom: carrier-composition predictor gbar_v for EVERY panel variant ----------
    # weight = founder baseline freq p0_f, restricted to analyzable founders with finite gamma
    per = {}   # chrom -> dict(pos,ref_len,alt_len,dlen,isins,maf, gbar_bio1,gbar_bio18, insload contribution)
    ins_load = np.zeros(founders.size); del_load = np.zeros(founders.size); tot_carry = np.zeros(founders.size)
    for ch in ("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"):
        cl = ch.lower(); base = f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        assert list(meta["founders"].astype("U6")) == list(founders), f"{ch} founder order mismatch"
        pos = meta["pos"].astype(np.int64); rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsr(); vc = sp.load_npz(f"{base}.var_called.npz").tocsr()
        n_alt = np.asarray(vp.sum(0)).ravel().astype(np.float64)
        n_cal = np.asarray(vc.sum(0)).ravel().astype(np.float64)
        vpc = vp.tocsc()
        dlen = np.abs(al - rl); isins = al > rl
        # founder-composition predictor per column (weighted mean gamma over analyzable carriers)
        gb = {}
        for k in ("bio1", "bio18"):
            w = np.where(analyz & np.isfinite(gamma[k]), p0f, 0.0)          # founder weights
            wg = w * np.nan_to_num(gamma[k])
            num = np.asarray((wg[None, :] @ vpc)).ravel()                  # sum_f w*gamma over carriers
            den = np.asarray((w[None, :] @ vpc)).ravel()                   # sum_f w over carriers
            gb[k] = np.divide(num, den, out=np.full(num.shape, np.nan), where=den > 0)
        maf = np.minimum(n_alt, n_cal - n_alt) / np.maximum(n_cal, 1)
        per[ch] = dict(pos=pos, rl=rl, al=al, dlen=dlen, isins=isins, maf=maf,
                       n_alt=n_alt, n_cal=n_cal, gb1=gb["bio1"], gb18=gb["bio18"])
        # founder SV insertion / deletion load (among common SVs, the analysis set)
        common = (n_alt >= MIN_MAC) & (n_alt <= vp.shape[0] - MIN_MAC) & (n_cal >= 0.9 * vp.shape[0])
        sv = common & (dlen > SV_BP)
        ins_cols = np.where(sv & isins)[0]; del_cols = np.where(sv & ~isins)[0]
        ins_load += np.asarray(vp[:, ins_cols].sum(1)).ravel()
        del_load += np.asarray(vp[:, del_cols].sum(1)).ravel()
        tot_carry += np.asarray(vp[:, np.where(common)[0]].sum(1)).ravel()
        svc = np.where(sv)[0]                        # capture SV-column keys+gbar only (for the T3 join)
        per[ch]["sv_keys"] = (ch, pos[svc], rl[svc], al[svc]); per[ch]["sv_g1"] = gb["bio1"][svc]; per[ch]["sv_g18"] = gb["bio18"][svc]
        print(f"[panel] {ch}: {vp.shape[1]:,} variants, common-SV ins={ins_cols.size:,} del={del_cols.size:,}", flush=True)

    # ================= T1: founder insertion-load vs climate response =================
    print("\n" + "=" * 74)
    print("T1  FOUNDER LEVEL: do founders carrying more insertions have more negative gamma?")
    a = analyz & np.isfinite(gamma["bio1"]) & (tot_carry > 0)
    ins_frac = ins_load / np.maximum(tot_carry, 1)      # insertion load normalized by variants carried
    for k in ("bio1", "bio18"):
        gk = gamma[k]
        r_cnt = stats.spearmanr(ins_load[a], gk[a]).correlation
        r_frac = stats.spearmanr(ins_frac[a], gk[a]).correlation
        r_del = stats.spearmanr(del_load[a], gk[a]).correlation
        print(f"  {k}: Spearman(gamma, insertion-load count)={r_cnt:+.3f}  "
              f"(ins-FRACTION)={r_frac:+.3f}  (deletion-load)={r_del:+.3f}   [n={int(a.sum())} founders]")
    print("  (negative r = founders with more insertions decline more as climate warms/dries;")
    print("   caveat: 231 founders are ~10-25 clades, so this is a FEW effective points.)")

    # ====== build a panel (chrom,pos,rl,al) -> mean gbar map over SV columns ONLY (for the T3 join) ======
    def panel_key(ch, pos, rl, al):
        return np.char.add(np.char.add(np.char.add(np.char.add(np.char.add(
            np.array([ch] * len(pos)), ":"), pos.astype(str)), ":"), rl.astype(str)),
            np.char.add(":", al.astype(str)))
    keymap_g1, keymap_g18 = {}, {}
    from collections import defaultdict
    acc1, acc18, accn = defaultdict(float), defaultdict(float), defaultdict(int)
    for ch, d in per.items():
        c, pos_s, rl_s, al_s = d["sv_keys"]
        ks = panel_key(c, pos_s, rl_s, al_s)
        for kk, g1, g18 in zip(ks, d["sv_g1"], d["sv_g18"]):
            acc1[kk] += (g1 if g1 == g1 else 0.0); acc18[kk] += (g18 if g18 == g18 else 0.0); accn[kk] += 1
    for kk in accn:
        keymap_g1[kk] = acc1[kk] / accn[kk]; keymap_g18[kk] = acc18[kk] / accn[kk]

    # ============ T2: does founder composition reproduce the SV-ins-vs-SNP / ins-vs-del shift? ============
    # panel-native SV set (ins/del) and SNP set, gbar computed above, matched on founder MAF
    print("\n" + "=" * 74)
    print("T2  Does founder-composition gbar (NO per-variant temporal data) reproduce the signal?")
    gb_all, maf_all, isins_all, issv_all, issnp_all = [], [], [], [], []
    for ch, d in per.items():
        common = (d["n_alt"] >= MIN_MAC) & (d["n_alt"] <= 231 - MIN_MAC) & (d["n_cal"] >= 0.9 * 231)
        sv = common & (d["dlen"] > SV_BP); snp = common & (d["rl"] == 1) & (d["al"] == 1)
        for k in ("bio1", "bio18"):
            pass
        gb_all.append(np.where(sv | snp, d["gb1"], np.nan))  # placeholder, we store per-axis below
    # store per-axis arrays
    store = {}
    for k, gbkey in (("bio1", "gb1"), ("bio18", "gb18")):
        GB, MAF, ISINS, ISSV, ISSNP = [], [], [], [], []
        for ch, d in per.items():
            common = (d["n_alt"] >= MIN_MAC) & (d["n_alt"] <= 231 - MIN_MAC) & (d["n_cal"] >= 0.9 * 231)
            sv = common & (d["dlen"] > SV_BP); snp = common & (d["rl"] == 1) & (d["al"] == 1)
            keep = (sv | snp) & np.isfinite(d[gbkey])
            GB.append(d[gbkey][keep]); MAF.append(d["maf"][keep]); ISINS.append(d["isins"][keep])
            ISSV.append(sv[keep]); ISSNP.append(snp[keep])
        GB = np.concatenate(GB); MAF = np.concatenate(MAF); ISINS = np.concatenate(ISINS)
        ISSV = np.concatenate(ISSV); ISSNP = np.concatenate(ISSNP)
        gb_sv = GB[ISSV]; gb_snp = GB[ISSNP]; maf_sv = MAF[ISSV]; maf_snp = MAF[ISSNP]
        ins = ISINS[ISSV]
        gb_ins, gb_del = gb_sv[ins], gb_sv[~ins]; maf_ins, maf_del = maf_sv[ins], maf_sv[~ins]
        m_snp_ins = match_to(maf_ins, gb_snp, maf_snp)
        print(f"\n  --- gbar predicted from d gamma / d {k} (founder composition only) ---")
        print(f"    median gbar: matchedSNP={np.median(m_snp_ins):+.5f}  SV-ins={np.median(gb_ins):+.5f}  "
              f"SV-del={np.median(gb_del):+.5f}   (n ins={ins.sum():,} del={(~ins).sum():,})")
        print(f"    KS SV-ins vs matched-SNP  p={stats.ks_2samp(gb_ins, m_snp_ins).pvalue:.2e}  "
              f"(frac gbar<0: ins={np.mean(gb_ins<0):.3f} mSNP={np.mean(m_snp_ins<0):.3f})")
        m_del_match = match_to(maf_del, gb_snp, maf_snp)
        print(f"    KS SV-del vs matched-SNP  p={stats.ks_2samp(gb_del, m_del_match).pvalue:.2e}  "
              f"(insertion vs deletion asymmetry present in founder-composition alone?)")
        store[f"gb_ins_{k}"] = gb_ins.astype(np.float32); store[f"gb_del_{k}"] = gb_del.astype(np.float32)

    # ============ T3: does gbar predict the ACTUAL temporal beta? align via callqual recipe ============
    print("\n" + "=" * 74)
    print("T3  Does founder-composition gbar predict the actual per-variant temporal beta?")
    b = np.load(f"{GEA}/sv_adaptive/s_climate_slope.npz")
    idx_non = np.load(f"{STORE}/index_nonsnp.npz")
    ch_non = idx_non["chrom"].astype("U5"); pos_non = idx_non["pos"].astype(np.int64)
    rl_n = idx_non["ref_len"].astype(np.int64); al_n = idx_non["alt_len"].astype(np.int64)
    dlen_n = np.abs(al_n - rl_n)
    common_non = lib.founder_panel_keep(ch_non, pos_non, min_mac=MIN_MAC)
    p0_non = np.load(f"{STORE}/p0_nonsnp.npy").astype(np.float64)
    non_m = (dlen_n >= 1) & common_non & (p0_non >= MIN_P0) & (p0_non <= 1 - MIN_P0)
    cols_non = np.where(non_m)[0]; sv_k = dlen_n[cols_non] > SV_BP
    ch_sv = ch_non[cols_non][sv_k]; pos_sv = pos_non[cols_non][sv_k]
    rl_sv = rl_n[cols_non][sv_k]; al_sv = al_n[cols_non][sv_k]
    assert ch_sv.size == b["p0_sv"].size, f"{ch_sv.size} vs {b['p0_sv'].size}"
    # per-row chrom (ch_sv already 'ChrN', U5) -> build keys directly (panel_key wants a scalar chrom)
    keys_sv = np.char.add(np.char.add(np.char.add(np.char.add(np.char.add(
        ch_sv.astype("U5"), ":"), pos_sv.astype(str)), ":"), rl_sv.astype(str)),
        np.char.add(":", al_sv.astype(str)))
    gsv1 = np.array([keymap_g1.get(kk, np.nan) for kk in keys_sv])
    gsv18 = np.array([keymap_g18.get(kk, np.nan) for kk in keys_sv])
    hit = np.isfinite(gsv1)
    print(f"    SV panel-join hit rate: {hit.mean():.3f}  (n={hit.sum():,}/{hit.size:,})")
    isdel_sv = b["sv_isdel"].astype(bool)
    for k, gsv in (("bio1", gsv1), ("bio18", gsv18)):
        beta = b[f"beta_sv_{k}"]; m = hit & np.isfinite(beta)
        rho = stats.spearmanr(gsv[m], beta[m]).correlation
        # conditional: residual beta after linear fit on gbar; does ins/del still differ?
        A = np.c_[np.ones(m.sum()), gsv[m]]; coef, *_ = np.linalg.lstsq(A, beta[m], rcond=None)
        resid = beta[m] - A @ coef
        di = isdel_sv[m]
        ks_raw = stats.ks_2samp(beta[m][~di], beta[m][di]).pvalue
        ks_res = stats.ks_2samp(resid[~di], resid[di]).pvalue
        print(f"    {k}: Spearman(gbar, beta)={rho:+.3f}   "
              f"ins-vs-del KS on RAW beta p={ks_raw:.1e}  ->  on beta|gbar RESIDUAL p={ks_res:.1e}   "
              f"(residual >> raw => asymmetry explained by founder composition)")

    np.savez_compressed(f"{GEA}/sv_adaptive/founder_load_test.npz",
                        gamma_bio1=gamma["bio1"], gamma_bio18=gamma["bio18"], founders=founders,
                        ins_load=ins_load, del_load=del_load, tot_carry=tot_carry, analyz=analyz, p0f=p0f,
                        gsv_bio1=gsv1, gsv_bio18=gsv18, beta_sv_bio1=b["beta_sv_bio1"],
                        beta_sv_bio18=b["beta_sv_bio18"], sv_isdel=b["sv_isdel"], **store)
    print("\n[wrote] founder_load_test.npz")


if __name__ == "__main__":
    main()
