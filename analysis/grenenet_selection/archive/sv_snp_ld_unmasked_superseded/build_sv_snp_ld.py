#!/usr/bin/env python
"""How much are our SVs tagged by SNPs? Per-SV max r2 with nearby SNPs.

Addresses the SNP-redundancy critique (Kang 2023: ~17% of Arabidopsis SVs tagged
at r2>0.6; Yan 2021: 54% of adaptive human SVs in strong SNP LD). For each true
SV (>50 bp) we compute the maximum genotype-LD r2 with any SNP within +/-WINDOW bp,
across the 231 founders. Two panels (run separately):

  panel      : SNPs from our own merged_231 panel (V_pa; same founders as the SVs)
  shortread  : SNPs from the GrENE-net short-read panel (greneNet_final_v1.1)

r2 = squared Pearson correlation of the two founder vectors (SV presence vs SNP
alt dosage), over the 231 shared founders. Monomorphic variants are skipped.

Output (--out): sv_snp_ld_<panel>_<chrom>.npz with per-SV chrom/pos/sv_size,
best_r2, best_snp_pos, n_snp_window.

This file computes the SV side from V_pa always; the SNP genotype matrix comes
from V_pa (panel) or from a prepared short-read genotype npz (shortread).
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

VPA = "panel/arch3/{cl}/var_pa_231_arch3_{cl}.var_pa.npz"
VMETA = "panel/arch3/{cl}/var_pa_231_arch3_{cl}.meta.npz"


def _load_vpa(chrom):
    cl = chrom.lower()
    z = np.load(VPA.format(cl=cl), allow_pickle=True)
    G = sp.csr_matrix((z["data"], z["indices"], z["indptr"]),
                      shape=tuple(z["shape"])).toarray().astype(np.float32)  # [231 x nrec]
    m = np.load(VMETA.format(cl=cl), allow_pickle=True)
    return G, m["pos"].astype(np.int64), m["ref_len"], m["alt_len"], \
        [str(f) for f in m["founders"]]


def _standardize(G):
    """[variants x founders] -> z-scored rows; returns (Z, ok_mask). Monomorphic -> 0."""
    mu = G.mean(1, keepdims=True); sd = G.std(1, keepdims=True)
    ok = (sd[:, 0] > 0)
    Z = np.zeros_like(G); Z[ok] = (G[ok] - mu[ok]) / sd[ok]
    return Z, ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--panel", choices=["panel", "shortread"], default="panel")
    ap.add_argument("--window", type=int, default=50000, help="+/- bp around SV")
    ap.add_argument("--sv-min-bp", type=int, default=50)
    ap.add_argument("--shortread-geno", default=None,
                    help="npz with snp_pos + geno [nSNP x 231] + founders (shortread mode)")
    ap.add_argument("--out", default=f"{lib.GEA}/sv_snp_ld")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    G, pos, ref_len, alt_len, founders = _load_vpa(args.chrom)   # G: [231 x nrec]
    n_founders = len(founders)
    sv = (np.abs(alt_len - ref_len) > args.sv_min_bp)
    sv_idx = np.where(sv)[0]
    sv_pos = pos[sv_idx]
    Gsv = G[:, sv_idx].T                                          # [nSV x 231]

    if args.panel == "panel":
        snpm = (ref_len == 1) & (alt_len == 1)
        snp_pos = pos[snpm]
        Gsnp = G[:, snpm].T                                      # [nSNP x 231]
        snp_founders = founders
    else:
        z = np.load(args.shortread_geno, allow_pickle=True)
        sr_founders = [str(f) for f in z["founders"]]
        # align short-read founders to our V_pa founder order
        idx = {f: i for i, f in enumerate(sr_founders)}
        keep = [f for f in founders if f in idx]
        col = [idx[f] for f in keep]
        Gsnp = z["geno"][:, col].astype(np.float32)             # [nSNP x len(keep)]
        snp_pos = z["snp_pos"].astype(np.int64)
        # subset SV founders to the same shared set
        svcol = [founders.index(f) for f in keep]
        Gsv = Gsv[:, svcol]
        n_founders = len(keep)
        print(f"shared founders: {len(keep)}/{len(founders)}", flush=True)

    # standardize once
    Zsv, sv_ok = _standardize(Gsv)
    Zsnp, snp_ok = _standardize(Gsnp)
    Zsnp = Zsnp[snp_ok]; snp_pos_ok = snp_pos[snp_ok]
    order = np.argsort(snp_pos_ok); snp_pos_ok = snp_pos_ok[order]; Zsnp = Zsnp[order]
    print(f"{args.chrom} {args.panel}: {len(sv_idx):,} SVs vs {Zsnp.shape[0]:,} "
          f"segregating SNPs (window +/-{args.window}, {n_founders} founders)", flush=True)

    best_r2 = np.full(len(sv_idx), np.nan, np.float32)
    best_pos = np.zeros(len(sv_idx), np.int64)
    n_win = np.zeros(len(sv_idx), np.int32)
    W = args.window
    lo = np.searchsorted(snp_pos_ok, sv_pos - W, "left")
    hi = np.searchsorted(snp_pos_ok, sv_pos + W, "right")
    for i in range(len(sv_idx)):
        if not sv_ok[i] or hi[i] <= lo[i]:
            continue
        block = Zsnp[lo[i]:hi[i]]                                 # [k x F]
        r2 = (block @ Zsv[i] / n_founders) ** 2                   # [k]
        n_win[i] = len(r2)
        j = int(np.argmax(r2))
        best_r2[i] = r2[j]; best_pos[i] = snp_pos_ok[lo[i] + j]
        if (i + 1) % 20000 == 0:
            print(f"  {i+1}/{len(sv_idx)}", flush=True)

    np.savez(f"{args.out}/sv_snp_ld_{args.panel}_{args.chrom}.npz",
             chrom=np.array([args.chrom] * len(sv_idx)), pos=sv_pos,
             sv_size=np.abs(alt_len - ref_len)[sv_idx], best_r2=best_r2,
             best_snp_pos=best_pos, n_snp_window=n_win)
    fin = np.isfinite(best_r2)
    print(f"{args.chrom} {args.panel}: SVs with a segregating SNP in window: "
          f"{fin.sum():,}/{len(sv_idx):,}")
    for t in (0.2, 0.5, 0.8, 0.95):
        print(f"  best_r2 > {t}: {(best_r2[fin] > t).mean()*100:.1f}%")
    print(f"-> {args.out}/sv_snp_ld_{args.panel}_{args.chrom}.npz")


if __name__ == "__main__":
    main()
