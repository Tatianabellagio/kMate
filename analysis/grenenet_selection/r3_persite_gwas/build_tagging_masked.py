#!/usr/bin/env python
"""build_tagging_masked.py -- properly missing-masked SV/indel-SNP tagging r2.

Fixes two bugs found while reviewing `build_sv_snp_ld.py`'s output in
sv_snp_ld_tagging.ipynb:

  1. Missingness was silently coded as REF. That script's `_load_vpa` only
     loads var_pa (ALT-carrier), never var_called (the called mask) -- so a
     missing genotype and a true REF call were indistinguishable. Panel SV
     records average ~60% missing (see panel_stats_arch3.ipynb SS4), and
     missingness is strongly cohort-correlated (long-read vs short-read), so
     this could inflate r2 for markers that are jointly missing in the same
     founders.
  2. No MAC floor. 53.5% of SVs are singletons (MAC=1) -- a singleton anchor
     can register a spurious r2=1 "tag" against any other marker private to
     the same one founder, a small-N coincidence rather than real LD.

Confirmed empirically that the GrENE-Net short-read VCF has ZERO missing
genotypes (sampled ~2M GT calls, no `./.`), so only the panel (arch3) side
needs call-mask handling; the anchor SV/indel side always needs it (its own
missingness), in both --panel modes.

Vectorized masked-correlation trick: the arch3 panel is haploid (0/1) for
every record, so the anchor's own variance always simplifies to
mean*(1-mean); the candidate/tag side uses the general (non-simplified) sum
of squares so it also works for the short-read side's 0/1/2 diploid dosage.
All four sufficient statistics per candidate window are matrix-vector
products (Sc @ v), not a per-SNP python loop -- this keeps the masked
version roughly as fast as the original unmasked one.

Usage (per chrom x anchor-class x panel-mode; run all 5 chroms x 2 classes
x 2 modes = 20 calls, ~shared loads if run per-chrom):

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY analysis/grenenet_selection/r3_persite_gwas/build_tagging_masked.py --chrom Chr1 \
      --anchor-class sv --panel-mode panel --mac-floor 2
"""
from __future__ import annotations
import argparse, os, sys, time
import numpy as np

PROJ = "/global/scratch/users/tbellg/kmate"
VPA   = f"{PROJ}/panel/arch3/{{cl}}/var_pa_231_arch3_{{cl}}.var_pa.npz"
VC    = f"{PROJ}/panel/arch3/{{cl}}/var_pa_231_arch3_{{cl}}.var_called.npz"
VMETA = f"{PROJ}/panel/arch3/{{cl}}/var_pa_231_arch3_{{cl}}.meta.npz"
SR_GENO = f"{PROJ}/analysis/grenenet_selection/sv_snp_ld/shortread_geno_{{Cl}}.npz"
OUTDIR = f"{PROJ}/analysis/grenenet_selection/sv_snp_ld_v2"


def _load_panel(chrom):
    cl = chrom.lower()
    zp = np.load(VPA.format(cl=cl), allow_pickle=False)
    zc = np.load(VC.format(cl=cl), allow_pickle=False)
    import scipy.sparse as sp
    G = sp.csr_matrix((zp["data"], zp["indices"], zp["indptr"]),
                       shape=tuple(zp["shape"])).toarray().astype(np.float32)   # [231 x N] dose (0/1)
    C = sp.csr_matrix((zc["data"], zc["indices"], zc["indptr"]),
                       shape=tuple(zc["shape"])).toarray().astype(np.float32)   # [231 x N] called (0/1)
    m = np.load(VMETA.format(cl=cl), allow_pickle=True)
    pos = m["pos"].astype(np.int64)
    rl = m["ref_len"].astype(np.int64); al = m["alt_len"].astype(np.int64)
    founders = [str(x) for x in m["founders"]]
    return G, C, pos, rl, al, founders


def classify(rl, al):
    ldiff = np.abs(al - rl)
    snp = (rl == 1) & (al == 1)
    sv = ldiff > 50
    indel = (~snp) & (~sv)
    cls = np.where(snp, "SNP", np.where(sv, "SV", "indel"))
    size = np.where(snp, 0, ldiff)
    return cls, size


def masked_r2_block(a_dose, a_call, S, Sc, min_pair):
    """a_dose,a_call: (F,). S,Sc: (K,F). Returns r2 (K,) with -1 sentinel for
    insufficient joint overlap, and n (K,) the joint-called count."""
    n = Sc @ a_call
    sum_a = Sc @ a_dose
    sum_s = S @ a_call
    sum_ss = (S * S) @ a_call
    sum_as = S @ a_dose
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_a = sum_a / n
        mean_s = sum_s / n
        var_a = mean_a * (1.0 - mean_a)          # anchor is always haploid 0/1
        var_s = sum_ss / n - mean_s ** 2          # general (handles 0/1/2 too)
        cov = sum_as / n - mean_a * mean_s
        denom = var_a * var_s
        r2 = np.where(denom > 1e-12, (cov ** 2) / denom, 0.0)
    r2 = np.where(n >= min_pair, r2, -1.0)
    return r2, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True, help="e.g. Chr1")
    ap.add_argument("--anchor-class", choices=["sv", "indel"], required=True)
    ap.add_argument("--panel-mode", choices=["panel", "shortread"], required=True)
    ap.add_argument("--window-bp", type=int, default=50_000)
    ap.add_argument("--mac-floor", type=int, default=2,
                     help="drop anchors/candidate-SNPs with MAC below this (default 2 = drop singletons)")
    ap.add_argument("--min-pair-frac", type=float, default=0.5,
                     help="min fraction of founders that must be jointly called to trust an r2")
    ap.add_argument("--out-dir", default=OUTDIR)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    t0 = time.time()

    G, C, pos, rl, al, founders = _load_panel(args.chrom)
    NF = len(founders)
    cls, size = classify(rl, al)
    print(f"[load] {args.chrom}: {len(pos):,} panel records ({NF} founders) in {time.time()-t0:.0f}s", file=sys.stderr)

    ac = G.sum(0); an = C.sum(0)
    maf = np.where(an > 0, np.minimum(ac, an - ac) / np.maximum(an, 1), 0.0)
    mac = np.minimum(ac, an - ac)
    mac_ok = mac >= args.mac_floor

    anchor_mask = (cls == ("SV" if args.anchor_class == "sv" else "indel")) & mac_ok
    anchor_idx = np.where(anchor_mask)[0]
    anchor_pos = pos[anchor_idx]
    print(f"[anchors] {args.anchor_class}: {mac_ok[cls == ('SV' if args.anchor_class=='sv' else 'indel')].sum():,} "
          f"of {(cls == ('SV' if args.anchor_class=='sv' else 'indel')).sum():,} pass MAC>={args.mac_floor} "
          f"-> {len(anchor_idx):,} anchors", file=sys.stderr)

    if args.panel_mode == "panel":
        snp_mask = (cls == "SNP") & mac_ok
        snp_idx = np.where(snp_mask)[0]
        snp_pos = pos[snp_idx]
        Sfull = G[:, snp_idx].T          # [nSNP x F]
        Scfull = C[:, snp_idx].T
        founder_order = founders
        n_shared = NF
    else:
        z = np.load(SR_GENO.format(Cl=args.chrom), allow_pickle=True)
        sr_founders = [str(x) for x in z["founders"]]
        idx_map = {f: i for i, f in enumerate(sr_founders)}
        shared = [f for f in founders if f in idx_map]
        if len(shared) < NF:
            print(f"  [warn] only {len(shared)}/{NF} founders shared with short-read VCF", file=sys.stderr)
        sr_col = [idx_map[f] for f in shared]
        panel_col = [founders.index(f) for f in shared]
        snp_pos = z["snp_pos"].astype(np.int64)
        Sfull = z["geno"][:, sr_col].astype(np.float32)     # [nSNP x n_shared], 0/1/2, fully called
        Scfull = np.ones_like(Sfull)
        # MAC floor on shortread SNPs (their own AC/AN over the shared founders)
        sr_ac = Sfull.sum(1) / 2.0  # rough allele count on a 0/1/2 scale -> /2 for allele units; fine for a floor
        sr_an = np.full(Sfull.shape[0], Sfull.shape[1], dtype=np.float64)
        sr_mac = np.minimum(sr_ac, sr_an - sr_ac)
        keep_snp = sr_mac >= args.mac_floor
        snp_pos, Sfull, Scfull = snp_pos[keep_snp], Sfull[keep_snp], Scfull[keep_snp]
        # restrict panel G/C to the shared founder columns for the anchor side
        G = G[panel_col]; C = C[panel_col]
        founder_order = shared
        n_shared = len(shared)
    order = np.argsort(snp_pos)
    snp_pos, Sfull, Scfull = snp_pos[order], Sfull[order], Scfull[order]
    print(f"[candidates] {args.panel_mode}: {len(snp_pos):,} SNPs pass MAC>={args.mac_floor}", file=sys.stderr)

    min_pair = int(args.min_pair_frac * n_shared)
    best_r2 = np.full(len(anchor_idx), np.nan, np.float32)
    best_pos = np.full(len(anchor_idx), -1, np.int64)
    n_win = np.zeros(len(anchor_idx), np.int32)
    W = args.window_bp
    lo_arr = np.searchsorted(snp_pos, anchor_pos - W, "left")
    hi_arr = np.searchsorted(snp_pos, anchor_pos + W, "right")

    if args.panel_mode == "panel":
        A_dose_all = G[:, anchor_idx]      # [F x nAnchor]
        A_call_all = C[:, anchor_idx]
    else:
        A_dose_all = G[:, anchor_idx]      # already restricted to shared founder rows above
        A_call_all = C[:, anchor_idx]

    t1 = time.time()
    for i in range(len(anchor_idx)):
        lo, hi = lo_arr[i], hi_arr[i]
        if hi <= lo:
            continue
        S = Sfull[lo:hi]; Sc = Scfull[lo:hi]
        a_dose = A_dose_all[:, i]; a_call = A_call_all[:, i]
        r2, n = masked_r2_block(a_dose, a_call, S, Sc, min_pair)
        valid = r2 >= 0
        n_win[i] = int(valid.sum())
        if not valid.any():
            continue
        j = np.argmax(np.where(valid, r2, -1))
        best_r2[i] = r2[j]
        best_pos[i] = snp_pos[lo:hi][j]
        if (i + 1) % 20000 == 0:
            print(f"  {args.anchor_class}/{args.panel_mode} {args.chrom}: {i+1:,}/{len(anchor_idx):,} "
                  f"({time.time()-t1:.0f}s)", file=sys.stderr)

    out = f"{args.out_dir}/tagging_{args.anchor_class}_{args.panel_mode}_{args.chrom}.npz"
    np.savez(out, chrom=np.array([args.chrom] * len(anchor_idx)), pos=anchor_pos,
             size=size[anchor_idx], ac=ac[anchor_idx], an=an[anchor_idx], maf=maf[anchor_idx],
             best_r2=best_r2, best_snp_pos=best_pos, n_snp_window=n_win)
    fin = np.isfinite(best_r2)
    print(f"[done] {args.anchor_class}/{args.panel_mode} {args.chrom}: {len(anchor_idx):,} anchors, "
          f"{fin.sum():,} with >=1 candidate; total {time.time()-t0:.0f}s -> {out}")
    for t in (0.2, 0.5, 0.8, 0.95):
        print(f"  r2>{t}: {100*(best_r2[fin] > t).mean():.1f}%")


if __name__ == "__main__":
    main()
