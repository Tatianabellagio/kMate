#!/usr/bin/env python
"""Spatial (block-rotation) null for the SV-on-selected-haplotype enrichment.

Addresses the audit's HIGH-1: the mac-matched cluster-resampling null treats hap-clusters as
independent, but the selected set is spatially clustered (per-block LD in p_joint) and sv_r2 is
spatially autocorrelated (SVs cluster in blocks), so the cluster-null under-disperses and p is
anti-conservative. The rotation null keeps the selected set's EXACT spatial pattern and slides
it along the genome against the FIXED sv_r2 landscape -- the standard co-localization null that
respects both spatial structures. Also reports how many DISTINCT BLOCKS drive the signal, and
re-tests every axis/MAF cell so multiplicity is explicit.

Env: kmate (MEMB_TAG=clq90). Writes sv_hap_rotationnull.csv.
"""
import os, sys
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
from founder_genotype import build_genotype

CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]; PANEL = "panel/arch3"
MAC_MIN = 3; nF = 231; NPERM = 10000
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
        na = np.asarray(vp.sum(0)).ravel(); ncl = np.asarray(vc.sum(0)).ravel()
        svj = np.where((dl > 50) & (na >= sv_mac) & (na <= nF - sv_mac) & (ncl / nF >= 0.9))[0]
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


def rotation_p(sv_r2_ord, sel_ord, rng, nperm=NPERM):
    """slide the selected set's spatial pattern along the genome; null = mean sv_r2 at shifts."""
    M = len(sv_r2_ord)
    obs = sv_r2_ord[sel_ord].mean()
    shifts = rng.integers(1, M, size=nperm)
    picks = (sel_ord[None, :] + shifts[:, None]) % M           # (nperm, n_sel)
    nul = sv_r2_ord[picks].mean(1)
    fold = obs / max(np.median(nul), 1e-12)
    p = (1 + (nul >= obs).sum()) / (nperm + 1)
    return obs, fold, p, np.median(nul)


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
    p_joint = raw["p_joint"]; p_clim = 2 * stats.norm.sf(np.abs(raw["z_clim"])); mac = raw["mac"].astype(int)
    M = len(p_joint)
    HG = pd.read_csv(HAPGEA)
    j = regp.merge(HG[["chrom", "unit_start", "unit_end", "cluster", "s_mean"]],
                   left_on=["chrom", "start", "end", "cluster"],
                   right_on=["chrom", "unit_start", "unit_end", "cluster"], how="left")
    s_abs = np.abs(j.s_mean.to_numpy())

    # genomic order for rotation
    chrom_code = regp.chrom.map({c: i for i, c in enumerate(CHROMS)}).to_numpy()
    gorder = np.lexsort((regp.start.to_numpy(), chrom_code))
    inv = np.empty(M, int); inv[gorder] = np.arange(M)          # cluster idx -> genomic rank
    blk_ord = regp.block.to_numpy()

    axes = {"JOINT": np.argsort(p_joint), "CLIMATE": np.argsort(p_clim), "TEMPORAL": np.argsort(-s_abs)}
    n1 = int(round(0.01 * M))
    rng = np.random.default_rng(0)
    rows = []
    print(f"{M:,} clusters | rotation null, {NPERM} shifts | top-1% = {n1}")
    print(f"{'axis':>9} {'MAC':>4} | {'obs':>7} {'rot_fold':>8} {'rot_p':>8} | "
          f"{'n_tag(r2>=.5)':>13} {'distinct_blocks':>15}")
    for svmac in (2, 6, 12, 24, 46):
        sv_r2 = sv_r2_for_mac(Gp, regp, svmac)
        sv_r2_ord = sv_r2[gorder]
        for axis in ("JOINT", "CLIMATE", "TEMPORAL"):
            sel = axes[axis][:n1]
            sel_ord = inv[sel]                                   # positions in genomic order
            obs, fold, p, med = rotation_p(sv_r2_ord, sel_ord, rng)
            # concentration: how many distinct blocks among selected TAGGING events
            tagsel = sel[sv_r2[sel] >= 0.5]
            ntag = len(tagsel); ndb = len(np.unique(blk_ord[tagsel]))
            print(f"{axis:>9} {svmac:>4} | {obs:>7.4f} x{fold:>6.2f} {p:>8.4f} | "
                  f"{ntag:>13} {ndb:>15}")
            rows.append(dict(axis=axis, sv_mac=svmac, obs=round(obs, 4), rot_fold=round(fold, 3),
                             rot_p=round(p, 4), n_tag_r2_05=ntag, distinct_blocks=ndb))
    df = pd.DataFrame(rows)
    df.to_csv(f"{lib.GEA}/sv_adaptive/sv_hap_rotationnull.csv", index=False)
    # multiplicity: BH across all cells
    df["bh_q"] = lib.bh(df.rot_p.to_numpy())
    nsig = int((df.rot_p < 0.05).sum()); nbh = int((df.bh_q < 0.1).sum())
    print(f"\ncells: {len(df)} | rot_p<0.05: {nsig} | BH q<0.1: {nbh}")
    print("surviving BH q<0.1:")
    print(df[df.bh_q < 0.1][["axis", "sv_mac", "rot_fold", "rot_p", "bh_q", "distinct_blocks"]].to_string(index=False))
    print(f"[wrote] {lib.GEA}/sv_adaptive/sv_hap_rotationnull.csv")


if __name__ == "__main__":
    main()
