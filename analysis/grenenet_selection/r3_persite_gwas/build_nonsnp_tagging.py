#!/usr/bin/env python
"""SNP-tagging r2 for every non-SNP marker that feeds K_nonsnp (indel+SV, MAC>=5, call>=90%).

For each such marker, the max founder-genotype r2 against any panel SNP within +/-50kb (all
segregating panel SNPs, no MAC floor on the SNP side -- gives SNPs their best shot at tagging).
Feeds the "SNP-untagged non-SNP GRM" test (VAREXP_SELECTION_HANDOFF.md next step): build a GRM
only from markers SNPs genuinely cannot see (best_r2 < 0.2), then ask whether that residual layer
explains selection-trait variance beyond K_snp.

Uses the SAME class definition + MAC/call filter as build_class_grms.py so the untagged marker
set is an exact subset of what feeds K_nonsnp -- carries the var_pa COLUMN INDEX (not just pos)
so downstream code can re-standardize the exact same markers with no position-matching ambiguity
(build_sv_snp_ld.py's prior SV-only run matched by position, fine there since it was a standalone
diagnostic; here we need exact GRM-column identity).

Output -> analysis/grenenet_selection/r3_persite_gwas/results/varexp/nonsnp_tagging_{chrom}.npz
  col_idx, pos, cls(1=indel,2=sv), best_r2, best_snp_pos, n_snp_window
Env: kmate. Run on a compute node (Lustre I/O + O(n_markers) windowed matmuls).
"""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r3_persite_gwas/results/varexp"
N = 231
MIN_MAC = 5
CALL_MIN = 0.9
WINDOW = 50_000


def _load_chrom(cl):
    base = f"{lib.PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
    meta = np.load(f"{base}.meta.npz", allow_pickle=True)
    vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
    vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
    n_alt = np.asarray(vp.sum(0)).ravel().astype(float)
    n_cal = np.asarray(vc.sum(0)).ravel().astype(float)
    return meta["pos"].astype(np.int64), meta["ref_len"].astype(int), meta["alt_len"].astype(int), \
        vp, vc, n_alt, n_cal


def _classes(rl, al):
    dlen = np.abs(al - rl)
    return np.where((rl == 1) & (al == 1), 0, np.where(dlen > 50, 2, 1))


def _dense_calledimputed(vp, vc, n_alt, n_cal, cols):
    g = vp[:, cols].toarray().astype(np.float32)
    called = vc[:, cols].toarray().astype(bool)
    p = np.where(n_cal[cols] > 0, n_alt[cols] / np.maximum(n_cal[cols], 1), 0.0)
    miss = ~called
    if miss.any():
        g[miss] = np.take(p, np.where(miss)[1])
    return g, p


def _standardize_cols(g, p):
    """g: [n x cols] -> standardized COLUMNS (variants), for a later variants-as-rows transpose."""
    sd = np.sqrt(np.maximum(p * (1 - p), 1e-12))
    ok = sd > 1e-6
    Z = np.zeros_like(g)
    Z[:, ok] = (g[:, ok] - p[ok]) / sd[ok]
    return Z.T, ok  # [cols x N]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True, choices=["chr1", "chr2", "chr3", "chr4", "chr5"])
    ap.add_argument("--window", type=int, default=WINDOW)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()

    pos, rl, al, vp, vc, n_alt, n_cal = _load_chrom(args.chrom)
    mac = np.minimum(n_alt, N - n_alt)
    keep = (mac >= MIN_MAC) & (n_cal >= CALL_MIN * N)
    vclass = _classes(rl, al)

    nonsnp_idx = np.where(keep & (vclass != 0))[0]
    snp_idx_all = np.where(vclass == 0)[0]           # SNP side: NOT MAC-floored (max tagging chance)
    print(f"{args.chrom}: {len(nonsnp_idx):,} non-SNP markers (indel="
          f"{(vclass[nonsnp_idx]==1).sum():,} sv={(vclass[nonsnp_idx]==2).sum():,}), "
          f"{len(snp_idx_all):,} candidate SNPs", flush=True)

    g_ns, p_ns = _dense_calledimputed(vp, vc, n_alt, n_cal, nonsnp_idx)
    Zns, ns_ok = _standardize_cols(g_ns, p_ns)        # [n_nonsnp x N]
    del g_ns

    g_sp, p_sp = _dense_calledimputed(vp, vc, n_alt, n_cal, snp_idx_all)
    Zsnp, sp_ok = _standardize_cols(g_sp, p_sp)       # [n_snp x N]
    del g_sp
    Zsnp = Zsnp[sp_ok]
    snp_pos_ok = pos[snp_idx_all][sp_ok]
    order = np.argsort(snp_pos_ok)
    snp_pos_ok = snp_pos_ok[order]
    Zsnp = Zsnp[order]

    ns_pos = pos[nonsnp_idx]
    best_r2 = np.full(len(nonsnp_idx), np.nan, np.float32)
    best_pos = np.zeros(len(nonsnp_idx), np.int64)
    n_win = np.zeros(len(nonsnp_idx), np.int32)
    W = args.window
    lo = np.searchsorted(snp_pos_ok, ns_pos - W, "left")
    hi = np.searchsorted(snp_pos_ok, ns_pos + W, "right")
    for i in range(len(nonsnp_idx)):
        if not ns_ok[i] or hi[i] <= lo[i]:
            continue
        block = Zsnp[lo[i]:hi[i]]                      # [k x N]
        r2 = (block @ Zns[i] / N) ** 2
        n_win[i] = len(r2)
        j = int(np.argmax(r2))
        best_r2[i] = r2[j]
        best_pos[i] = snp_pos_ok[lo[i] + j]
        if (i + 1) % 50_000 == 0:
            print(f"  {i+1:,}/{len(nonsnp_idx):,} ({time.time()-t0:.0f}s)", flush=True)

    out = f"{args.out}/nonsnp_tagging_{args.chrom}.npz"
    np.savez(out, col_idx=nonsnp_idx, pos=ns_pos, cls=vclass[nonsnp_idx],
             best_r2=best_r2, best_snp_pos=best_pos, n_snp_window=n_win)
    fin = np.isfinite(best_r2)
    print(f"{args.chrom}: markers with a segregating SNP in window: {fin.sum():,}/{len(nonsnp_idx):,}")
    for t in (0.2, 0.5, 0.8, 0.95):
        print(f"  best_r2 > {t}: {(best_r2[fin] > t).mean()*100:.1f}%")
    print(f"untagged (best_r2 < 0.2, incl. no-SNP-in-window): "
          f"{(~(best_r2 > 0.2)).sum():,}/{len(nonsnp_idx):,}")
    print(f"-> {out}  ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
