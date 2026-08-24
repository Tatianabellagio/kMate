#!/usr/bin/env python
"""Per-variant-class founder GRMs for the SNP-vs-non-SNP variance partition.

Builds standardized genomic-relationship matrices (231 x 231) from the arch3 founder carrier
matrix var_pa, one per variant CLASS, so we can partition the ecotype-selection trait's variance
into what SNPs vs non-SNP variation (indels + SVs) explain:

  K_snp     SNPs                (ref_len==1 & alt_len==1)
  K_indel   small non-SNP       (|alt_len-ref_len| <= 50)
  K_sv      structural          (|alt_len-ref_len| >  50)
  K_nonsnp  indel + SV pooled   (the phase-1-vs-kMate contrast)
  K_all     everything
  K_snp_matched  SNPs subsampled to the non-SNP COUNT and MAC spectrum (fairness control:
                 is any SNP advantage just "more markers / different MAF"?)

Each K = Z Z' / M with Z = (g - 2p)/sqrt(2p(1-p)) style standardization on the founder allele
frequency p (uncalled cells imputed to p); markers filtered to MAC>=MIN_MAC & call>=CALL_MIN.
GRMs are built on all 231 founders (relatedness is a founder property); the modeling step subsets
to the analyzable set. Two streaming passes over the 5 per-chrom var_pa npz (pass 1: class GRMs +
MAC spectra; pass 2: MAC-matched SNP GRM). Heavy Lustre I/O -> run via sbatch, NOT the login node.

Output -> analysis/grenenet_selection/varexp/class_grms.npz
  founders, K_snp, K_indel, K_sv, K_nonsnp, K_all, K_snp_matched, n_markers{dict}, MIN_MAC, CALL_MIN
Env: kmate.
"""
from __future__ import annotations
import os, sys, time, json
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/varexp"
CHROMS = ["chr1", "chr2", "chr3", "chr4", "chr5"]
N = 231
MIN_MAC = 5           # MAF ~2% floor for GRM markers
CALL_MIN = 0.9
BLK = 200_000
NBIN = 20             # MAC bins for the matched-SNP control
SEED = 0


def _load_chrom(cl):
    base = f"{lib.PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
    meta = np.load(f"{base}.meta.npz", allow_pickle=True)
    vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
    vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
    n_alt = np.asarray(vp.sum(0)).ravel().astype(float)
    n_cal = np.asarray(vc.sum(0)).ravel().astype(float)
    return (meta["ref_len"].astype(int), meta["alt_len"].astype(int), vp, vc, n_alt, n_cal)


def _classes(rl, al):
    """0=SNP, 1=indel(|dlen|<=50), 2=SV(|dlen|>50)."""
    dlen = np.abs(al - rl)
    return np.where((rl == 1) & (al == 1), 0, np.where(dlen > 50, 2, 1))


def _ZZt(vp, vc, n_alt, n_cal, cols):
    """Accumulate Z Z' over `cols` (uncalled -> column freq p). Returns (231x231, len(cols))."""
    K = np.zeros((N, N))
    for s in range(0, len(cols), BLK):
        c = cols[s:s + BLK]
        g = vp[:, c].toarray().astype(np.float64)
        called = vc[:, c].toarray().astype(bool)
        p = np.where(n_cal[c] > 0, n_alt[c] / np.maximum(n_cal[c], 1), 0.0)
        miss = ~called
        if miss.any():
            g[miss] = np.take(p, np.where(miss)[1])
        sd = np.sqrt(np.maximum(p * (1 - p), 1e-6))   # 0/1 founder presence (inbred): haploid coding
        Z = (g - p) / sd
        K += Z @ Z.T
    return K, len(cols)


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    founders = np.load(f"{lib.PROJ}/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz",
                       allow_pickle=True)["founders"].astype(str)

    # ---- pass 1: per-class GRMs + non-SNP / SNP MAC spectra ----
    Ksum = {k: np.zeros((N, N)) for k in ("snp", "indel", "sv")}
    Mc = {k: 0 for k in ("snp", "indel", "sv")}
    nonsnp_mac, snp_mac_bychrom = [], {}
    print("pass 1: class GRMs", flush=True)
    for cl in CHROMS:
        rl, al, vp, vc, n_alt, n_cal = _load_chrom(cl)
        mac = np.minimum(n_alt, N - n_alt)
        keep = (mac >= MIN_MAC) & (n_cal >= CALL_MIN * N)
        vclass = _classes(rl, al)
        for name, cid in (("snp", 0), ("indel", 1), ("sv", 2)):
            cols = np.where(keep & (vclass == cid))[0]
            K, m = _ZZt(vp, vc, n_alt, n_cal, cols)
            Ksum[name] += K; Mc[name] += m
        nonsnp_mac.append(mac[keep & (vclass != 0)])
        snp_mac_bychrom[cl] = (np.where(keep & (vclass == 0))[0], mac)   # keep for pass 2
        print(f"  {cl}: snp+={Mc['snp']:,} indel+={Mc['indel']:,} sv+={Mc['sv']:,} "
              f"({time.time()-t0:.0f}s)", flush=True)

    K_snp = Ksum["snp"] / Mc["snp"]
    K_indel = Ksum["indel"] / Mc["indel"]
    K_sv = Ksum["sv"] / Mc["sv"]
    K_nonsnp = (Ksum["indel"] + Ksum["sv"]) / (Mc["indel"] + Mc["sv"])
    K_all = (Ksum["snp"] + Ksum["indel"] + Ksum["sv"]) / (Mc["snp"] + Mc["indel"] + Mc["sv"])

    # ---- pass 2: MAC-matched SNP GRM (match non-SNP count & MAC spectrum) ----
    print("pass 2: MAC-matched SNP GRM", flush=True)
    nonsnp_mac = np.concatenate(nonsnp_mac)
    edges = np.unique(np.quantile(nonsnp_mac, np.linspace(0, 1, NBIN + 1)))
    tgt, _ = np.histogram(nonsnp_mac, bins=edges)                  # non-SNP count per MAC bin
    # total SNPs per bin across chroms
    snp_bin_tot = np.zeros(len(tgt))
    for cl in CHROMS:
        cols, mac = snp_mac_bychrom[cl]
        snp_bin_tot += np.histogram(mac[cols], bins=edges)[0]
    p_keep = np.divide(tgt, np.maximum(snp_bin_tot, 1), where=snp_bin_tot > 0)
    p_keep = np.clip(p_keep, 0, 1)
    rng = np.random.default_rng(SEED)
    Km, Mm = np.zeros((N, N)), 0
    for cl in CHROMS:
        cols, mac = snp_mac_bychrom[cl]
        b = np.clip(np.digitize(mac[cols], edges) - 1, 0, len(tgt) - 1)
        sub = cols[rng.random(len(cols)) < p_keep[b]]
        rl, al, vp, vc, n_alt, n_cal = _load_chrom(cl)
        K, m = _ZZt(vp, vc, n_alt, n_cal, sub); Km += K; Mm += m
        print(f"  {cl}: matched SNPs += {Mm:,} ({time.time()-t0:.0f}s)", flush=True)
    K_snp_matched = Km / Mm

    nmark = {"snp": Mc["snp"], "indel": Mc["indel"], "sv": Mc["sv"],
             "nonsnp": Mc["indel"] + Mc["sv"], "all": sum(Mc.values()),
             "snp_matched": Mm}
    np.savez(f"{OUT}/class_grms.npz", founders=founders,
             K_snp=K_snp, K_indel=K_indel, K_sv=K_sv, K_nonsnp=K_nonsnp, K_all=K_all,
             K_snp_matched=K_snp_matched, MIN_MAC=MIN_MAC, CALL_MIN=CALL_MIN,
             n_markers=json.dumps(nmark))
    print("\nn_markers:", json.dumps(nmark))
    print("GRM diag means:", {k: round(float(np.diag(v).mean()), 2) for k, v in
                              [("snp", K_snp), ("nonsnp", K_nonsnp), ("sv", K_sv),
                               ("snp_matched", K_snp_matched)]})
    print(f"corr(K_snp, K_nonsnp) offdiag = "
          f"{np.corrcoef(K_snp[np.triu_indices(N,1)], K_nonsnp[np.triu_indices(N,1)])[0,1]:.3f}")
    print(f"wrote {OUT}/class_grms.npz  ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
