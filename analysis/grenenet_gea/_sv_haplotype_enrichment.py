#!/usr/bin/env python
"""CORRECT unit test: is the SIGNAL-CARRYING HAPLOTYPE (hap-cluster) SV-tagged?

The published sv_enrichment.py works at the BLOCK level: rank a block by min-p_joint over its
hap-clusters, then count SVs anywhere in the block. That flags "block contains a selected
haplotype" AND "block contains an SV (on ANY of its ~5 haplotypes)" -- it never checks the SV
is on the SELECTED haplotype. This aligns the unit to the founder-GWAS's actual test unit
(the hap-cluster): rank HAP-CLUSTERS by p_joint, and ask whether a selected hap-cluster is
more likely to TAG a common SV (r²(SV, cluster) >= thr, SV positioned in the cluster's block)
than a founder-frequency(mac)-matched non-selected hap-cluster.

If the block-level ×2.07 was real at the haplotype level, selected hap-clusters should be
SV-tagged more than matched controls. If it collapses, the block-level number was co-occurrence.
Env: kmate (MEMB_TAG=clq90).
"""
import os, sys
import numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from founder_genotype import build_genotype

CHROMS = ["Chr1","Chr2","Chr3","Chr4","Chr5"]; PANEL = "panel/arch3"
MAC_MIN = 3; nF = 231
SV_MAC = int(os.environ.get("SV_MAC", 12))   # SV founder floor: 12=common(MAF5%), 2=no-singletons


def r2_cols(a, B):
    a = a - a.mean(); Bc = B - B.mean(0)
    num = (a[:, None] * Bc).sum(0) ** 2; den = (a @ a) * (Bc ** 2).sum(0)
    return np.divide(num, den, out=np.zeros_like(num), where=den > 0)


def main():
    # founder-GWAS markers = hap-clusters (replicate the poly filter so Gp aligns to raw Z)
    G, founders, reg = build_genotype()                 # MEMB_TAG=clq90
    reg["unit"] = reg.chrom + ":" + reg.start.astype(str) + "-" + reg.end.astype(str)
    cnt = G.sum(0); blk = reg.block.to_numpy()
    order = np.lexsort((-cnt, blk)); sbk = blk[order]
    is_ref = np.zeros(len(cnt), bool); is_ref[order[np.concatenate([[True], sbk[1:] != sbk[:-1]])]] = True
    poly = (~is_ref) & (cnt >= MAC_MIN) & (cnt <= nF - MAC_MIN)
    Gp = G[:, poly].astype(float); regp = reg[poly].reset_index(drop=True)
    raw = lib.multisite_gwas_raw("clq90_pc1")
    assert list(regp.unit.values) == list(raw["unit"]), "Gp/raw marker alignment failed"
    p_joint = raw["p_joint"]; mac = raw["mac"].astype(int)   # cluster founder count
    M = len(p_joint)
    print(f"{M:,} hap-cluster markers (the founder-GWAS unit) | {regp.unit.nunique():,} blocks")

    # SV-tagging: for each hap-cluster, max r² to a common SV positioned in its block
    sv_r2 = np.zeros(M)
    unit_rows = {u: g.index.to_numpy() for u, g in regp.groupby("unit")}
    for ch in CHROMS:
        cl = ch.lower()
        meta = np.load(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.meta.npz", allow_pickle=True)
        pos = meta["pos"].astype(np.int64); dl = np.abs(meta["alt_len"].astype(np.int64) - meta["ref_len"].astype(np.int64))
        vp = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_pa.npz")
        vc = sp.load_npz(f"{PANEL}/{cl}/var_pa_231_arch3_{cl}.var_called.npz")
        na = np.asarray(vp.sum(0)).ravel(); nc = np.asarray(vc.sum(0)).ravel()
        svj = np.where((dl > 50) & (na >= SV_MAC) & (na <= nF - SV_MAC) & (nc / nF >= 0.9))[0]
        # block bounds per unit on this chrom
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
    print(f"hap-clusters that tag a common SV: r²>=0.5 {int((sv_r2>=0.5).sum()):,} | "
          f"r²>=0.9 {int((sv_r2>=0.9).sum()):,}  (of {M:,})")

    # matched enrichment on mac (cluster founder-frequency) bins
    edges = np.unique(np.quantile(mac, np.linspace(0, 1, 21)).astype(int))
    binid = np.clip(np.digitize(mac, edges[1:-1]), 0, len(edges) - 2)
    order_p = np.argsort(p_joint)
    print("\nCORRECT haplotype-level enrichment: P(selected hap-cluster tags an SV) vs mac-matched null")
    print(f"{'top':>6} {'n':>6} | {'tag r²>=0.5':>22} | {'tag r²>=0.9':>22}")
    rng = np.random.default_rng(0)
    csv_rows = []
    for frac in (0.005, 0.01, 0.02):
        n = int(round(frac * M)); sel = order_p[:n]
        out = []
        for thr in (0.5, 0.9):
            tag = (sv_r2 >= thr).astype(float)
            obs = tag[sel].mean()
            # size(mac)-matched null
            nullm = np.empty(2000)
            binmem = {b: np.where(binid == b)[0] for b in np.unique(binid[sel])}
            selbin = binid[sel]
            for k in range(2000):
                pick = [rng.choice(binmem[b]) for b in selbin]
                nullm[k] = tag[pick].mean()
            fold = obs / max(np.median(nullm), 1e-9)
            p = (1 + (nullm >= obs).sum()) / 2001
            out.append(f"x{fold:.2f} p={p:.3f} ({obs:.3f}v{np.median(nullm):.3f})")
            csv_rows.append(dict(sv_mac=SV_MAC, top=f"{frac:.1%}", n_selected=n, r2_thr=thr,
                                 obs_tag_rate=round(obs, 4), null_tag_rate=round(float(np.median(nullm)), 4),
                                 fold=round(fold, 3), p_perm=round(p, 4)))
        print(f"{frac:>6.1%} {n:>6} | {out[0]:>22} | {out[1]:>22}")
    outcsv = f"{lib.GEA}/sv_adaptive/sv_haplotype_enrichment_mac{SV_MAC}.csv"
    pd.DataFrame(csv_rows).to_csv(outcsv, index=False)
    print(f"[wrote] {outcsv}")

    # continuous: max r²-to-SV of selected clusters vs matched
    print("\ncontinuous: median max-r²(cluster,SV) among top clusters vs all clusters:")
    for frac in (0.005, 0.01, 0.02):
        n = int(round(frac * M)); sel = order_p[:n]
        print(f"  top {frac:.1%}: selected {np.median(sv_r2[sel]):.3f} | genome-wide {np.median(sv_r2):.3f} | "
              f"mean selected {sv_r2[sel].mean():.3f} vs {sv_r2.mean():.3f}")


if __name__ == "__main__":
    main()
