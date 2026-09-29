#!/usr/bin/env python
"""Is the frequency-dependent SV enrichment real, or a frequency-alignment artifact?

Hypothesis (user): r²(SV,haplotype) is maximized at equal frequency, and the selection
statistics have frequency biases (JOINT->common haplotypes, TEMPORAL->rare via noisier
|s_mean|). When the SV MAF floor and the selection ranking point at the same-frequency
haplotypes, coarse mac-matching leaves a spurious enrichment.

Test: (1) report each axis's frequency bias (median cluster mac of top-1% vs all);
(2) re-run the key significant cells with EXACT-mac matched nulls (null drawn only from
haplotypes of the IDENTICAL founder count). If the signal survives exact matching it is real;
if it collapses, the whole frequency-dependence was the artifact. Env: kmate (MEMB_TAG=clq90).
"""
import os, sys
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "r3_persite_gwas"))
import lib
from founder_genotype import build_genotype

CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
PANEL = "panel/arch3"
MAC_MIN = 3
nF = 231
NPERM = 2000
HAPGEA = "results/grenenet_gea/hapfreq_clq90/pipelineB_varlen/hap_gea.csv"


def r2_cols(a, B):
    a = a - a.mean(); Bc = B - B.mean(0)
    num = (a[:, None] * Bc).sum(0) ** 2; den = (a @ a) * (Bc ** 2).sum(0)
    return np.divide(num, den, out=np.zeros_like(num), where=den > 0)


def sv_r2_for_mac(Gp, regp, sv_mac):
    M = Gp.shape[1]; sv_r2 = np.zeros(M)
    for ch in CHROMS:
        cl = ch.lower()
        meta = np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64)
        dl = np.abs(meta["alt_len"].astype(np.int64) - meta["ref_len"].astype(np.int64))
        vp = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz")
        vc = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz")
        na = np.asarray(vp.sum(0)).ravel(); nc = np.asarray(vc.sum(0)).ravel()
        svj = np.where((dl > 50) & (na >= sv_mac) & (na <= nF - sv_mac) & (nc / nF >= 0.9))[0]
        rc = regp[regp.chrom == ch]
        for u, g in rc.groupby("unit"):
            s0, e0 = g.start.iloc[0], g.end.iloc[0]
            inb = svj[(pos[svj] >= s0) & (pos[svj] <= e0)]
            if len(inb) == 0:
                continue
            cols = g.index.to_numpy(); Bc = Gp[:, cols]
            for j in inb:
                a = vp[:, j].toarray().ravel().astype(float)
                sv_r2[cols] = np.maximum(sv_r2[cols], r2_cols(a, Bc))
    return sv_r2


def exact_mac_null(tagvec, sel, mac, rng, nperm=NPERM):
    """null mean(tagvec) drawing each selected cluster's control from IDENTICAL-mac clusters."""
    by = {m: np.where(mac == m)[0] for m in np.unique(mac[sel])}
    tot = np.zeros(nperm)
    for i in sel:
        pool = by[mac[i]]
        tot += tagvec[pool[rng.integers(0, len(pool), nperm)]]
    return tot / len(sel)


def main():
    os.chdir("/global/scratch/users/tbellg/kmate")
    G, founders, reg = build_genotype()
    reg["unit"] = reg.chrom + ":" + reg.start.astype(str) + "-" + reg.end.astype(str)
    reg["cluster"] = reg.groupby("block").cumcount()
    cnt = G.sum(0); blk = reg.block.to_numpy()
    order = np.lexsort((-cnt, blk)); sbk = blk[order]
    is_ref = np.zeros(len(cnt), bool)
    is_ref[order[np.concatenate([[True], sbk[1:] != sbk[:-1]])]] = True
    poly = (~is_ref) & (cnt >= MAC_MIN) & (cnt <= nF - MAC_MIN)
    Gp = G[:, poly].astype(float); regp = reg[poly].reset_index(drop=True)
    raw = lib.multisite_gwas_raw("clq90_pc1")
    assert list(regp.unit.values) == list(raw["unit"])
    p_joint = raw["p_joint"]; p_clim = 2 * stats.norm.sf(np.abs(raw["z_clim"]))
    mac = raw["mac"].astype(int); M = len(p_joint)
    HG = pd.read_csv(HAPGEA)
    j = regp.merge(HG[["chrom", "unit_start", "unit_end", "cluster", "s_mean"]],
                   left_on=["chrom", "start", "end", "cluster"],
                   right_on=["chrom", "unit_start", "unit_end", "cluster"], how="left")
    s_abs = np.abs(j.s_mean.to_numpy())
    axes = {"JOINT": np.argsort(p_joint), "CLIMATE": np.argsort(p_clim), "TEMPORAL": np.argsort(-s_abs)}

    # (1) frequency bias of each axis
    print("=== frequency bias: cluster founder-count (mac) of top-1% vs all ===")
    print(f"{'axis':>9}  median_mac_top1%  median_mac_all  mean_top1%  mean_all")
    n1 = int(round(0.01 * M))
    for a, o in axes.items():
        sel = o[:n1]
        print(f"{a:>9}  {np.median(mac[sel]):>15.0f}  {np.median(mac):>13.0f}  "
              f"{mac[sel].mean():>9.1f}  {mac.mean():>7.1f}")

    # (2) exact-mac vs coarse-bin for the key cells
    edges = np.unique(np.quantile(mac, np.linspace(0, 1, 21)).astype(int))
    binid = np.clip(np.digitize(mac, edges[1:-1]), 0, len(edges) - 2)
    uidx = np.arange(M)

    def coarse_null(tagvec, sel, rng):
        selbin = binid[sel]; tot = np.zeros(NPERM)
        for b, c in zip(*np.unique(selbin, return_counts=True)):
            mem = uidx[binid == b]
            tot += tagvec[mem[rng.integers(0, len(mem), size=(NPERM, int(c)))]].sum(1)
        return tot / len(sel)

    cells = [("JOINT", 24), ("JOINT", 46), ("TEMPORAL", 2), ("TEMPORAL", 6), ("CLIMATE", 46)]
    print("\n=== continuous mean-r² enrichment: coarse-bin vs EXACT-mac null (top 1%) ===")
    print(f"{'axis':>9} {'MAC':>4} | {'obs':>7} | {'coarse fold/p':>18} | {'EXACT fold/p':>18}  verdict")
    rng = np.random.default_rng(0)
    rows = []
    macs_needed = sorted(set(m for _, m in cells))
    svr2 = {m: sv_r2_for_mac(Gp, regp, m) for m in macs_needed}
    for axis, svmac in cells:
        sel = axes[axis][:n1]; tag = svr2[svmac]
        obs = tag[sel].mean()
        cn = coarse_null(tag, sel, rng); en = exact_mac_null(tag, sel, mac, rng)
        cf, cp = obs / max(np.median(cn), 1e-12), (1 + (cn >= obs).sum()) / (NPERM + 1)
        ef, ep = obs / max(np.median(en), 1e-12), (1 + (en >= obs).sum()) / (NPERM + 1)
        verdict = "SURVIVES" if ep < 0.05 else "collapses -> artifact"
        print(f"{axis:>9} {svmac:>4} | {obs:>7.4f} | x{cf:>5.2f} p={cp:>6.4f} | "
              f"x{ef:>5.2f} p={ep:>6.4f}  {verdict}")
        rows.append(dict(axis=axis, sv_mac=svmac, obs=round(obs, 4), coarse_fold=round(cf, 3),
                         coarse_p=round(cp, 4), exact_fold=round(ef, 3), exact_p=round(ep, 4),
                         survives_exact=(ep < 0.05)))
    pd.DataFrame(rows).to_csv(f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive/sv_hap_freqrobust.csv", index=False)
    print(f"\n[wrote] {lib.GEA}/r1_sv_negative_selection/results/sv_adaptive/sv_hap_freqrobust.csv")


if __name__ == "__main__":
    main()
