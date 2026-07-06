"""Tier 2 at the CORRECT unit (user 2026-07-02): are the FITNESS-WINNING HAPLOTYPES SV-tagged?

The block test was wrong: a fitness block carrying an SV somewhere is meaningless if the SV is
not on the SELECTED haplotype. So test at the haplotype (hap-cluster) unit, exactly the unit the
selection signal is defined on ([[feedback-test-at-signal-unit]]).

Unit = clq90 hap-cluster (blocks_mcf90/hap_membership/chr{N}_hapmemb_clq90.npz: labels[unit,231]
= founder->cluster). A haplotype = a set of member founders.
  HAP FITNESS (per axis)  = sum_{f in cluster} founder Delta_h  = the cluster's frequency change
                            (the ecotype-sorting signal, kinship-relevant, from ecotype_fitness).
  SELECTED (winner)       = top-q genome-wide of |hap fitness|; zone-union = max over
                            {cold,mid,hot}x{rel,cen} (conditional neutrality), global separately.
  SV-TAGS-HAP             = max over common in-unit SVs of r2(SV carrier[231], cluster indicator
                            [231]).  i.e. is a common SV in LD with THIS haplotype across founders.
  NULL                    = EXACT founder-count matched (r2 is frequency-aligned, so each winner
                            hap's controls have the IDENTICAL member count). This is the control
                            the other agent's regional covariates do NOT provide, and it uses the
                            kinship-corrected garden-fitness axis instead of the JOINT founder-GWAS.

SV MAF floors swept (MAC 12/24/46 = MAF 5/10/20%) to reproduce the frequency dependence.

Output -> results/grenenet_gea/ecotype_fitness/gwas/hap_sv_enrichment.{csv,json}
Env: kmate.
"""
from __future__ import annotations
import os, sys, json
import numpy as np, pandas as pd, scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

G = f"{lib.GEA}/ecotype_fitness/gwas"
MEMB = f"{lib.GEA}/blocks_mcf90/hap_membership"
# HAPMEMB_TAG selects the founder->hap partition: "clq90" (default, all-class) or "clq90nosv"
# (SNP+indel-only clustering, self-tagging circularity control). SFX suffixes all outputs so the
# nosv run does not clobber the original table.
MEMB_TAG = os.environ.get("HAPMEMB_TAG", "clq90")
SFX = "" if MEMB_TAG == "clq90" else f"_{MEMB_TAG}"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
PANEL = "panel/arch3"
nF = 231
ZONE_AX = ["w_cold_rel", "w_mid_rel", "w_hot_rel", "w_cold_cen", "w_mid_cen", "w_hot_cen"]
GLOB_AX = ["w_global_rel", "w_global_cen"]
SV_MACS = [12, 24, 46]
TOP = [0.01, 0.005]
NPERM = 5000


def r2_vec(a, B):
    """r2 between vector a (n,) and each column of B (n,k), over founders."""
    a = a - a.mean(); Bc = B - B.mean(0)
    num = (a[:, None] * Bc).sum(0) ** 2
    den = (a @ a) * (Bc ** 2).sum(0)
    return np.divide(num, den, out=np.zeros_like(num), where=den > 0)


def main():
    # founder fitness aligned to membership founder order
    fit = pd.read_csv(f"{lib.GEA}/ecotype_fitness/ecotype_fitness.csv").set_index("founder")
    m0 = np.load(f"{MEMB}/chr1_hapmemb_{MEMB_TAG}.npz", allow_pickle=True)
    founders = m0["founders"].astype(str)
    fit = fit.loc[[int(f) if f.isdigit() else f for f in founders]] \
        if fit.index.map(str).isin(founders).any() else fit
    # robust align by string
    fidx = {str(k): i for i, k in enumerate(fit.index.astype(str))}
    order = [fidx[f] for f in founders]
    F = {ax: fit.iloc[order][ax].to_numpy(float) for ax in ZONE_AX + GLOB_AX}

    haps = []   # (chrom, start, end, ncount, tag_by_mac{mac:..}, fit_zone, fit_glob)
    for ch in CHROMS:
        cl = ch.lower()
        z = np.load(f"{MEMB}/{cl}_hapmemb_{MEMB_TAG}.npz", allow_pickle=True)
        labels = z["labels"]; us, ue = z["unit_start"], z["unit_end"]; uk = z["unit_k"]
        meta = np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64)
        dl = np.abs(meta["alt_len"].astype(np.int64) - meta["ref_len"].astype(np.int64))
        vp = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz").tocsc()
        vc = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz").tocsc()
        na = np.asarray(vp.sum(0)).ravel(); nc = np.asarray(vc.sum(0)).ravel()
        sv_ok = {mac: (dl > 50) & (na >= mac) & (na <= nF - mac) & (nc / nF >= 0.9)
                 for mac in SV_MACS}
        for u in range(len(labels)):
            lab = labels[u]
            s0, e0 = us[u], ue[u]
            in_unit = np.where((pos >= s0) & (pos <= e0))[0]
            # common SV carriers in this unit, per mac floor
            sv_by_mac = {mac: in_unit[sv_ok[mac][in_unit]] for mac in SV_MACS}
            for c in range(uk[u]):
                ind = (lab == c).astype(float)          # cluster indicator over 231 founders
                k = int(ind.sum())
                if k < 2 or k > nF - 2:
                    continue
                tag = {}
                for mac in SV_MACS:
                    js = sv_by_mac[mac]
                    if len(js) == 0:
                        tag[mac] = 0.0
                    else:
                        A = vp[:, js].toarray().astype(float)   # 231 x nsv
                        tag[mac] = float(r2_vec(ind, A).max())
                fz = max(abs(np.dot(ind, F[ax])) for ax in ZONE_AX)   # |sum member Δh| per axis
                fg = max(abs(np.dot(ind, F[ax])) for ax in GLOB_AX)
                haps.append((ch, int(s0), int(e0), k, tag[12], tag[24], tag[46], fz, fg))
        print(f"{ch}: {len(haps)} haps so far", flush=True)

    H = pd.DataFrame(haps, columns=["chrom", "start", "end", "k",
                                    "tag12", "tag24", "tag46", "fit_zone", "fit_glob"])
    H.to_csv(f"{G}/hap_sv_table{SFX}.csv", index=False)

    def exact_k_null(tag, kcol, sel_idx, rng, nperm=NPERM):
        obs = tag[sel_idx].mean()
        by_k = {}
        for kk in np.unique(kcol[sel_idx]):
            by_k[kk] = np.where(kcol == kk)[0]
        null = np.empty((len(sel_idx), nperm))
        for i, s in enumerate(sel_idx):
            pool = by_k[kcol[s]]
            null[i] = tag[pool][rng.integers(0, len(pool), nperm)]
        nm = null.mean(0)
        fold = obs / nm.mean() if nm.mean() > 0 else np.nan
        p = (1 + (nm >= obs).sum()) / (nperm + 1)
        return float(obs), float(nm.mean()), float(fold), float(p)

    rng = np.random.default_rng(0)
    kcol = H.k.to_numpy()
    rows = []
    for axis, fcol in (("zone_union", H.fit_zone.to_numpy()), ("global", H.fit_glob.to_numpy())):
        for q in TOP:
            thr = np.quantile(fcol, 1 - q)
            sel = np.where(fcol >= thr)[0]
            for mac, tname in ((12, "tag12"), (24, "tag24"), (46, "tag46")):
                tag = H[tname].to_numpy()
                obs, nm, fold, p = exact_k_null(tag, kcol, sel, rng)
                rows.append(dict(axis=axis, top_q=q, sv_mac=mac, n_sel=len(sel),
                                 obs_r2=round(obs, 4), null_r2=round(nm, 4),
                                 fold=round(fold, 3), p=round(p, 4)))
    out = pd.DataFrame(rows)
    out.to_csv(f"{G}/hap_sv_enrichment{SFX}.csv", index=False)
    json.dump(rows, open(f"{G}/hap_sv_enrichment{SFX}.json", "w"), indent=2)
    print(f"\nhaplotypes tested: {len(H):,}  (k in [2,{nF-2}])")
    print("\n=== winning haplotypes SV-tagged? (kinship-corrected fitness axis, exact-k null) ===")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
