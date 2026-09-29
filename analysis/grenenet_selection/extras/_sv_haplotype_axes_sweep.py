#!/usr/bin/env python
"""Correct-UNIT SV enrichment across ALL THREE selection axes + MAF sweep + power.

Closes the gap left by _sv_haplotype_enrichment.py, which only tested the FITNESS/JOINT axis
at the haplotype unit. Here every axis is scored at the founder-GWAS's actual unit -- the
hap-cluster -- and asked the SAME question: is a SELECTED hap-cluster more likely to TAG a
common SV (max r²(SV, cluster) >= thr, SV positioned in the cluster's block) than a
founder-frequency(mac)-matched non-selected hap-cluster?

Axes (all on the identical hap-cluster universe, mac-matched null):
  JOINT     -- founder GWAS any-site selection   (rank by p_joint asc)
  CLIMATE   -- founder GWAS climate contrast      (rank by p_clim asc, from z_clim)
  TEMPORAL  -- pool-seq generational selection    (rank by |s_mean| desc, from hap_gea.csv)

MAF sweep: SV founder floor SV_MAC in {2,6,12,24,46} (MAF ~0.9/2.6/5/10/20%) -- does the
conclusion move with the floor? (For JOINT we already know MAC 2 and 12 are both null.)

POWER: per cell, the minimum enrichment fold detectable at p<0.05 (null's 95th pct / median).
POSITIVE CONTROL: spike a KNOWN fold into the JOINT ranking and confirm the test recovers it
-- proves a real haplotype-level enrichment would not be missed.

Env: kmate (MEMB_TAG=clq90). Writes sv_adaptive/sv_haplotype_axes_sweep.csv + _power.csv +
_poscontrol.csv.
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
SV_MACS = [2, 6, 12, 24, 46]
FRACS = (0.005, 0.01, 0.02)
THRS = (0.5, 0.9)
NPERM = 2000
HAPGEA = "results/grenenet_gea/hapfreq_clq90/pipelineB_varlen/hap_gea.csv"


def r2_cols(a, B):
    a = a - a.mean()
    Bc = B - B.mean(0)
    num = (a[:, None] * Bc).sum(0) ** 2
    den = (a @ a) * (Bc ** 2).sum(0)
    return np.divide(num, den, out=np.zeros_like(num), where=den > 0)


def sv_r2_for_mac(Gp, regp, sv_mac):
    """max r²(SV, cluster) over common SVs (founder count in [sv_mac, nF-sv_mac]) in each
    cluster's block."""
    M = Gp.shape[1]
    sv_r2 = np.zeros(M)
    n_sv_used = 0
    for ch in CHROMS:
        cl = ch.lower()
        meta = np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64)
        dl = np.abs(meta["alt_len"].astype(np.int64) - meta["ref_len"].astype(np.int64))
        vp = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz")
        vc = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz")
        na = np.asarray(vp.sum(0)).ravel()
        nc = np.asarray(vc.sum(0)).ravel()
        svj = np.where((dl > 50) & (na >= sv_mac) & (na <= nF - sv_mac) & (nc / nF >= 0.9))[0]
        n_sv_used += len(svj)
        rc = regp[regp.chrom == ch]
        for u, g in rc.groupby("unit"):
            s0, e0 = g.start.iloc[0], g.end.iloc[0]
            inb = svj[(pos[svj] >= s0) & (pos[svj] <= e0)]
            if len(inb) == 0:
                continue
            cols = g.index.to_numpy()
            Bc = Gp[:, cols]
            for j in inb:
                a = vp[:, j].toarray().ravel().astype(float)
                rr = r2_cols(a, Bc)
                sv_r2[cols] = np.maximum(sv_r2[cols], rr)
    return sv_r2, n_sv_used


def matched_null_mean(tagvec, sel, binid, universe_idx, rng, nperm=NPERM):
    """null distribution of mean(tagvec) over mac-matched random picks of |sel| clusters.

    Vectorized: for each mac bin, draw (nperm x n_in_bin) picks at once and accumulate the
    tag sum per permutation. Reproducible (seeded rng), fast (no per-cluster Python loop)."""
    selbin = binid[sel]
    ubins, counts = np.unique(selbin, return_counts=True)
    tot = np.zeros(nperm)
    for b, c in zip(ubins, counts):
        mem = universe_idx[binid[universe_idx] == b]
        idx = rng.integers(0, len(mem), size=(nperm, int(c)))
        tot += tagvec[mem[idx]].sum(axis=1)
    return tot / len(sel)


def main():
    os.chdir("/global/scratch/users/tbellg/kmate")
    G, founders, reg = build_genotype()                       # MEMB_TAG=clq90
    reg["unit"] = reg.chrom + ":" + reg.start.astype(str) + "-" + reg.end.astype(str)
    reg["cluster"] = reg.groupby("block").cumcount()
    cnt = G.sum(0)
    blk = reg.block.to_numpy()
    order = np.lexsort((-cnt, blk))
    sbk = blk[order]
    is_ref = np.zeros(len(cnt), bool)
    is_ref[order[np.concatenate([[True], sbk[1:] != sbk[:-1]])]] = True
    poly = (~is_ref) & (cnt >= MAC_MIN) & (cnt <= nF - MAC_MIN)
    Gp = G[:, poly].astype(float)
    regp = reg[poly].reset_index(drop=True)

    raw = lib.multisite_gwas_raw("clq90_pc1")
    assert list(regp.unit.values) == list(raw["unit"]), "Gp/raw marker alignment failed"
    p_joint = raw["p_joint"]
    z_clim = raw["z_clim"]
    p_clim = 2 * stats.norm.sf(np.abs(z_clim))
    mac = raw["mac"].astype(int)
    M = len(p_joint)

    # temporal score per cluster, aligned by (chrom,start,end,cluster)
    HG = pd.read_csv(HAPGEA)
    j = regp.merge(HG[["chrom", "unit_start", "unit_end", "cluster", "n_founders", "s_mean"]],
                   left_on=["chrom", "start", "end", "cluster"],
                   right_on=["chrom", "unit_start", "unit_end", "cluster"], how="left")
    assert j.s_mean.notna().all(), "temporal join incomplete"
    assert (j.n_founders.to_numpy() == mac).all(), "temporal mac mismatch"
    s_abs = np.abs(j.s_mean.to_numpy())

    axes = {
        "JOINT":    np.argsort(p_joint),           # ascending p
        "CLIMATE":  np.argsort(p_clim),
        "TEMPORAL": np.argsort(-s_abs),            # descending |s_mean|
    }
    print(f"{M:,} hap-cluster markers | axes: JOINT, CLIMATE, TEMPORAL | MAF sweep {SV_MACS}")

    # mac bins shared across axes (same universe = all M)
    edges = np.unique(np.quantile(mac, np.linspace(0, 1, 21)).astype(int))
    binid = np.clip(np.digitize(mac, edges[1:-1]), 0, len(edges) - 2)
    universe_idx = np.arange(M)
    rng = np.random.default_rng(0)

    rows, prows = [], []
    for sv_mac in SV_MACS:
        sv_r2, n_sv = sv_r2_for_mac(Gp, regp, sv_mac)
        base05 = (sv_r2 >= 0.5).mean()
        base09 = (sv_r2 >= 0.9).mean()
        print(f"\nSV_MAC={sv_mac} (MAF~{sv_mac/(2*nF)*2*100:.1f}%): {n_sv:,} common SVs | "
              f"base tag-rate r²>=0.5 {base05:.4f}, r²>=0.9 {base09:.4f} | "
              f"mean r² {sv_r2.mean():.4f}")
        for axis, ordr in axes.items():
            for frac in FRACS:
                n = int(round(frac * M))
                sel = ordr[:n]
                # threshold metrics
                for thr in THRS:
                    tag = (sv_r2 >= thr).astype(float)
                    obs = tag[sel].mean()
                    nul = matched_null_mean(tag, sel, binid, universe_idx, rng)
                    med = max(np.median(nul), 1e-9)
                    fold = obs / med
                    p = (1 + (nul >= obs).sum()) / (NPERM + 1)
                    q95 = np.quantile(nul, 0.95)
                    min_fold = q95 / med                       # detectable fold at p~0.05
                    rows.append(dict(sv_mac=sv_mac, axis=axis, top=f"{frac:.1%}", n_sel=n,
                                     metric=f"tag_r2>={thr}", obs=round(obs, 4),
                                     null=round(med, 4), fold=round(fold, 3), p_perm=round(p, 4),
                                     min_detect_fold=round(min_fold, 3),
                                     n_events=int(round(obs * n))))
                # continuous mean r²
                obs_c = sv_r2[sel].mean()
                nulc = matched_null_mean(sv_r2, sel, binid, universe_idx, rng)
                medc = max(np.median(nulc), 1e-12)
                p_c = (1 + (nulc >= obs_c).sum()) / (NPERM + 1)
                rows.append(dict(sv_mac=sv_mac, axis=axis, top=f"{frac:.1%}", n_sel=n,
                                 metric="mean_r2", obs=round(obs_c, 5), null=round(medc, 5),
                                 fold=round(obs_c / medc, 3), p_perm=round(p_c, 4),
                                 min_detect_fold=round(np.quantile(nulc, 0.95) / medc, 3),
                                 n_events=n))
        # per-mac power summary line
        prows.append(dict(sv_mac=sv_mac, n_sv=n_sv, base_tag_r2_05=round(base05, 4),
                          base_tag_r2_09=round(base09, 4), mean_r2=round(float(sv_r2.mean()), 4)))

    df = pd.DataFrame(rows)
    out = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive/sv_haplotype_axes_sweep.csv"
    df.to_csv(out, index=False)
    pd.DataFrame(prows).to_csv(f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive/sv_haplotype_axes_power.csv", index=False)
    print(f"\n[wrote] {out}")

    # ---- positive control: spike a known fold into the JOINT top-1%, confirm recovery ----
    print("\n=== POSITIVE CONTROL (spike-in on JOINT top-1%, SV_MAC=12 base rate) ===")
    sv_r2_12, _ = sv_r2_for_mac(Gp, regp, 12)
    base = (sv_r2_12 >= 0.5).mean()
    ordr = axes["JOINT"]
    n = int(round(0.01 * M))
    sel = ordr[:n]
    sel_mask = np.zeros(M, bool); sel_mask[sel] = True
    pc_rows = []
    for target in (1.0, 1.5, 2.0, 3.0):
        rng2 = np.random.default_rng(123)
        # synthetic tags: non-selected at base rate, selected at target*base
        p_tag = np.where(sel_mask, np.clip(target * base, 0, 1), base)
        tag = (rng2.random(M) < p_tag).astype(float)
        obs = tag[sel].mean()
        nul = matched_null_mean(tag, sel, binid, universe_idx, rng2)
        med = max(np.median(nul), 1e-9)
        p = (1 + (nul >= obs).sum()) / (NPERM + 1)
        pc_rows.append(dict(spiked_fold=target, recovered_fold=round(obs / med, 3),
                            obs=round(obs, 4), null=round(med, 4), p_perm=round(p, 4),
                            n_events=int(round(obs * n)), detected=(p < 0.05)))
        print(f"  spiked x{target}: recovered x{obs/med:.2f}  p={p:.4f}  "
              f"({int(round(obs*n))} tagged of {n})  {'DETECTED' if p<0.05 else 'missed'}")
    pd.DataFrame(pc_rows).to_csv(f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive/sv_haplotype_axes_poscontrol.csv", index=False)
    print("[wrote] sv_haplotype_axes_poscontrol.csv")


if __name__ == "__main__":
    main()
