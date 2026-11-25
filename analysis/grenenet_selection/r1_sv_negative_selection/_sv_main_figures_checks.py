#!/usr/bin/env python
"""Robustness checks behind `sv_main_figures.ipynb` (events, truvari collapse defaults).

1. LOAD NULL (Fig. 2c, 3a). The per-garden Wilcoxon treats the 31 gardens as independent, but they
   share one founder pool. Null here: 1000 random sets of shared CDS SNPs with the SAME carrier-count x
   carrier-climate cells as the SV class (one SNP drawn per event, within its cell); statistic = mean
   over gardens of the log load change of the set minus that of the re-weighted matched SNPs (as the
   notebook). Empirical one-sided p (two-sided for deletions reported too).
2. FOUNDER ORIGIN (Fig. 2a-b). The burden LMM refit with founder origin bio1 (and bio1^2) as extra
   covariates, and Spearman burden-fitness within each origin tertile.
3. NO SINGLETONS (Fig. 1a, d). MAF spectra with minor-allele count >= 2: share MAF < 5% and KS vs synonymous SNPs.

Env: kmate. Writes results/sv_adaptive/main_fig_checks_{load,lmm,nosingleton}.csv.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
import importlib.util
_H = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location("lmm", os.path.join(_H, "_sv_burden_lmm.py"))
lm = importlib.util.module_from_spec(_s); _s.loader.exec_module(lm)

G = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
NDRAW = 1000
CALL_MIN = float(os.environ.get("SV_CALL_MIN", 0.9))   # founder call-rate floor; != 0.9 -> outputs tagged _cr{pct}
CRTAG = "" if CALL_MIN == 0.9 else f"_cr{int(round(CALL_MIN * 100))}"
CB = np.array([2, 3, 4, 5, 6, 8, 11, 16, 26, 51, 1000])


def main():
    rng = np.random.default_rng(0)
    # ---------- 1. load null
    E0 = pd.read_csv(f"{G}/sv_event_catalog.csv"); P0 = np.load(f"{G}/sv_event_catalog_p.npy")
    ok = ((E0.n_carriers >= 2) & (E0.n_called >= CALL_MIN * 231)).to_numpy(); E = E0[ok].reset_index(drop=True); PE = P0[ok].astype(float)
    GS = np.load(f"{G}/grenenet_sfs_event{CRTAG}.npz", allow_pickle=True); PV = GS["p"]; KD = GS["kind"]; CX = GS["ctx"]
    BQ = np.quantile(GS["b1"][KD == "snp"], [0.2, 0.4, 0.6, 0.8])
    CELL = np.digitize(GS["n"], CB) * 10 + np.digitize(GS["b1"], BQ); SH = (GS["nc"] > 0) & (GS["nw"] > 0)
    ms = np.where(SH & (KD == "snp") & (CX == 0))[0]; cellS = CELL[ms]; PS = PV[ms].astype(float)
    cn = pd.Series(cellS).value_counts(); bycell = {c: np.where(cellS == c)[0] for c in cn.index}
    cellE = np.digitize(E.n_carriers, CB) * 10 + np.digitize(E.carrier_bio1, BQ)

    def stat(Ls, cells):
        w = (pd.Series(cells).value_counts() / cn).reindex(cellS).fillna(0).to_numpy()
        Ln = (PS * w[:, None]).sum(0)
        return (np.log(Ls[1:] / Ls[0]) - np.log(Ln[1:] / Ln[0])).mean()
    C = E[(E.context == "CDS") & (E.origin == "shared")]
    CI_, CD_ = C[C.kind == "ins"], C[C.kind == "del"]
    classes = [("all CDS insertions", CI_.index), ("TE-derived", CI_.index[CI_.seq_class == "te_derived"]),
               ("other duplication", CI_.index[CI_.seq_class.isin(["dispersed_dup", "local_dup"])]),
               ("gene duplicate", CI_.index[CI_.seq_class == "gene_dup"]), ("novel", CI_.index[CI_.seq_class == "novel"]),
               ("all CDS deletions", CD_.index)]
    rows = []
    for lab, idx in classes:
        idx = np.asarray(idx); cells = cellE[idx]; kp = np.isin(cells, cn.index); idx, cells = idx[kp], cells[kp]
        obs = stat(PE[idx].sum(0), cells)
        null = np.array([stat(PS[np.array([rng.choice(bycell[c]) for c in cells])].sum(0), cells) for _ in range(NDRAW)])
        rows.append(dict(cls=lab, n_shared=idx.size, load=obs, null_mean=null.mean(), null_sd=null.std(),
                         p_lower=(np.sum(null <= obs) + 1) / (NDRAW + 1),
                         p_two_sided=(np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean())) + 1) / (NDRAW + 1)))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(f"{G}/main_fig_checks_load{CRTAG}.csv", index=False)

    # ---------- 2. founder origin in the burden LMM
    D = pd.read_csv(f"{G}/founder_burden_event{CRTAG}.csv", dtype={"founder": str}); K = np.load(f"{G}/founder_kinship.npy")
    Zf = np.load(lm.EF, allow_pickle=True)
    y = np.log((Zf["H0"][None, :] + Zf["DH_rel"] + lm.EPS) / (Zf["H0"][None, :] + lm.EPS)).mean(0); L = lm.LMM(K); z = lm.z
    rate = lambda c: (D[c] / D[f"called_{c}"]).to_numpy(float)
    t, r = rate("ins_shared_CDS_all"), rate("ins_shared_intergenic_all")
    base = [np.ones(len(D)), z(D.callrate_sv), z(D.callrate_snp), D.is_cactus.astype(float)]
    rows = []
    for lab, extra in (("as in the notebook", []), ("+ founder origin bio1", [z(D.bio1)]), ("+ bio1 + bio1^2", [z(D.bio1), z(D.bio1) ** 2])):
        b, se, _ = L.fit(y, np.column_stack(base + extra + [z(t), z(r)]))
        rows.append(dict(check=lab, beta=b[-2], se=se[-2], p=2 * stats.norm.sf(abs(b[-2] / se[-2])), n=len(D)))
    for g in ("cold", "temperate", "warm"):
        mk = (D.group == g).to_numpy(); rr = stats.spearmanr(t[mk], y[mk])
        rows.append(dict(check=f"Spearman within {g}-origin founders", beta=rr.statistic, se=np.nan, p=rr.pvalue, n=int(mk.sum())))
    pd.DataFrame(rows).to_csv(f"{G}/main_fig_checks_lmm{CRTAG}.csv", index=False); print(pd.DataFrame(rows).round(4).to_string(index=False))

    # ---------- 3. spectra without singletons
    CS = np.load(f"{G}/classic_sfs{CRTAG}.npz")
    def maf(n, minc):
        c = CS[f"{n}_c"].astype(float); k = CS[f"{n}_k"].astype(float); kp = (c >= minc) & (k - c >= minc); f = c[kp] / k[kp]
        return np.minimum(f, 1 - f)
    rows = []
    for minc in (1, 2):
        ref = maf("snp_syn", minc)
        for n in ("snp_syn", "snp_nonsyn", "ins_CDS", "ins_UTR", "ins_intron", "ins_TE", "ins_promoter", "ins_intergenic", "del_CDS"):
            v = maf(n, minc); ks = stats.ks_2samp(v, ref)
            rows.append(dict(min_minor_count=minc, cls=n, n=v.size, maf_lt5=(v < 0.05).mean(), KS_D=ks.statistic, KS_p=ks.pvalue))
    pd.DataFrame(rows).to_csv(f"{G}/main_fig_checks_nosingleton{CRTAG}.csv", index=False)
    print(pd.DataFrame(rows).round(4).to_string(index=False))
    print(f"[wrote] {G}/main_fig_checks_{{load,lmm,nosingleton}}{CRTAG}.csv")


if __name__ == "__main__":
    main()
