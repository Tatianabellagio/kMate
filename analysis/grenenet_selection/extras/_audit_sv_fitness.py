#!/usr/bin/env python
"""AUDIT the claim: common SVs (MAF>=10%) enriched on garden-fitness-selected haplotypes.

Independent re-derivation + sanity checks, deliberately NOT reusing the production r2/null code:
  1. ALIGNMENT   -- does Gp column i map to the same hap-cluster as raw p_joint[i]?
                    (mac from the genotype matrix must equal the GWAS npz mac, element-wise)
  2. DIRECTION   -- argsort(p_joint) top = largest chi2 = most selected?
  3. REIMPL      -- recompute sv_r2 with np.corrcoef (not r2_cols) and the enrichment with an
                    exact-mac resample written from scratch; must reproduce ~x1.9-2.0, p~0.01.
  4. SHUFFLE     -- permute p_joint across clusters -> enrichment fold must collapse to ~1.
Env: kmate (MEMB_TAG=clq90).
"""
import os, sys
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
from founder_genotype import build_genotype

CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]; PANEL = "panel/arch3"
MAC_MIN = 3; nF = 231; NPERM = 2000


def sv_r2_corrcoef(Gp, regp, sv_mac):
    """INDEPENDENT sv_r2: max r^2 via np.corrcoef (different code path from r2_cols)."""
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
            cols = g.index.to_numpy(); B = Gp[:, cols]
            for j in inb:
                a = vp[:, j].toarray().ravel().astype(float)
                if a.std() == 0:
                    continue
                for k, ci in enumerate(cols):
                    b = B[:, k]
                    if b.std() == 0:
                        continue
                    r = np.corrcoef(a, b)[0, 1]
                    sv_r2[ci] = max(sv_r2[ci], r * r)
    return sv_r2


def exact_resample(tag, sel, mac, rng, nperm=NPERM):
    """from-scratch exact-mac null; explicitly EXCLUDE each selected cluster from its own pool."""
    idx = np.arange(len(tag))
    tot = np.zeros(nperm)
    for i in sel:
        pool = idx[(mac == mac[i]) & (idx != i)]
        if len(pool) == 0:
            pool = idx[mac == mac[i]]
        tot += tag[pool[rng.integers(0, len(pool), nperm)]]
    return tot / len(sel)


def main():
    os.chdir("/global/scratch/users/tbellg/kmate")
    G, founders, reg = build_genotype()
    reg["unit"] = reg.chrom + ":" + reg.start.astype(str) + "-" + reg.end.astype(str)
    cnt = G.sum(0); blk = reg.block.to_numpy()
    order = np.lexsort((-cnt, blk)); sbk = blk[order]
    is_ref = np.zeros(len(cnt), bool)
    is_ref[order[np.concatenate([[True], sbk[1:] != sbk[:-1]])]] = True
    poly = (~is_ref) & (cnt >= MAC_MIN) & (cnt <= nF - MAC_MIN)
    Gp = G[:, poly].astype(float); regp = reg[poly].reset_index(drop=True)
    raw = lib.multisite_gwas_raw("clq90_pc1")
    p_joint = raw["p_joint"]; chi2 = raw["chi2_joint"]; mac_raw = raw["mac"].astype(int)
    M = len(p_joint)

    # 1. ALIGNMENT
    mac_G = Gp.sum(0).astype(int)
    agree = float((mac_G == mac_raw).mean())
    print(f"[1 ALIGNMENT] Gp mac == GWAS npz mac : {agree*100:.4f}% of {M:,} clusters "
          f"({'PASS' if agree==1.0 else 'FAIL -- p_joint mis-mapped to clusters!'})")
    if agree < 1.0:
        bad = np.where(mac_G != mac_raw)[0][:5]
        print("   first mismatches (idx, mac_G, mac_raw):", [(int(i), int(mac_G[i]), int(mac_raw[i])) for i in bad])

    # 2. DIRECTION
    o = np.argsort(p_joint)
    print(f"[2 DIRECTION] top-10 by argsort(p_joint): chi2 median {np.median(chi2[o[:10]]):.1f} vs "
          f"overall median {np.median(chi2):.1f} ; p range {p_joint[o[0]]:.2e}..{p_joint[o[9]]:.2e} "
          f"({'PASS: top=most selected' if np.median(chi2[o[:10]])>np.median(chi2) else 'FAIL: inverted'})")

    # 3. INDEPENDENT REIMPL (MAC 24 and 46)
    n1 = int(round(0.01 * M))
    for svmac in (24, 46):
        svr2 = sv_r2_corrcoef(Gp, regp, svmac)
        sel = o[:n1]; obs = svr2[sel].mean()
        rng = np.random.default_rng(7)
        nul = exact_resample(svr2, sel, mac_raw, rng)
        fold = obs / max(np.median(nul), 1e-12); p = (1 + (nul >= obs).sum()) / (NPERM + 1)
        print(f"[3 REIMPL MAC{svmac}] independent corrcoef sv_r2 + from-scratch exact-mac null: "
              f"obs {obs:.4f}  fold x{fold:.2f}  p={p:.4f}")
        # 4. SHUFFLE (same svmac)
        folds = []
        rs = np.random.default_rng(0)
        for _ in range(20):
            ps = rs.permutation(M)                      # shuffle cluster<->p_joint mapping
            sels = ps[:n1]
            nul2 = exact_resample(svr2, sels, mac_raw, rng)
            folds.append(svr2[sels].mean() / max(np.median(nul2), 1e-12))
        folds = np.array(folds)
        print(f"[4 SHUFFLE MAC{svmac}] 20 label-shuffles: mean fold {folds.mean():.2f} "
              f"(sd {folds.std():.2f}, max {folds.max():.2f}) "
              f"({'PASS: collapses to ~1' if folds.mean()<1.2 else 'FAIL: shuffled still enriched'})")


if __name__ == "__main__":
    main()
