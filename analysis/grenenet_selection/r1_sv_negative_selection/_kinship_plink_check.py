#!/usr/bin/env python
"""Cross-check the founder GRM: our NumPy construction vs plink2 --make-rel, on the SAME SNPs.

`_sv_burden_build.py` builds the kinship matrix in NumPy -- standardise each SNP, K = Z Z' / m --
rather than with a genotype toolkit. That is the GCTA/VanRaden-style GRM, but "we wrote it ourselves"
is a fair thing for a reviewer to poke at, so this reproduces the exact SNP set (same filters, same
seed, same chromosome order -- the RNG is used only for this subsample) and computes the same matrix
two ways:

  ours    Z = (G - mean) / sd  per SNP, K = Z Z' / m           (missing calls filled with the SNP mean)
  plink2  --make-rel square, i.e. the standard realised relationship matrix, on a VCF written from the
          same genotypes with the accessions coded homozygous (0/0, 1/1) as inbred lines

The two differ in one respect that matters for interpretation: we standardise by the OBSERVED SD of
each SNP, plink2 by sqrt(2p(1-p)) under Hardy-Weinberg. For inbred homozygous lines these differ by a
constant factor per SNP, so the matrices should agree up to scale -- the off-diagonal correlation is
the number to look at, not the absolute values.

env: kmate (+ plink2 from gwas_tools). Compute node, ~10 min, ~8 GB.
Writes results/sv_adaptive/kinship_check/{grm_ours.npy,plink2.rel,comparison.csv}.
"""
from __future__ import annotations
import os, sys, subprocess
import numpy as np, pandas as pd
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import lib                                                       # noqa: E402

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
CHK = f"{OUT}/kinship_check"
PLINK2 = os.path.expanduser("~/miniforge3/envs/gwas_tools/bin/plink2")
SEED = 0          # as in _sv_burden_build.py
FRAC = 0.02       # 2% SNP subsample
os.makedirs(CHK, exist_ok=True)


def collect_snps():
    """The kinship SNP set of _sv_burden_build.py: biallelic SNPs, called in >= 90% of founders,
    5 <= carriers <= N-5, then a 2% draw -- one RNG stream, chromosomes in lib.CHROMS order."""
    rng = np.random.default_rng(SEED)
    FO = pd.read_csv(f"{OUT}/sv_origin_partition_founders.csv", dtype={"founder": str})
    N = len(FO)
    blocks, ids = [], []
    for C in lib.CHROMS:
        base = f"{lib.PROJ}/panel/arch3/{C.lower()}/var_pa_231_arch3_{C.lower()}"
        meta = np.load(f"{base}.meta.npz", allow_pickle=True)
        rl = meta["ref_len"].astype(np.int64); al = meta["alt_len"].astype(np.int64)
        pos = meta["pos"].astype(np.int64)
        vp = sp.load_npz(f"{base}.var_pa.npz").tocsc()
        vc = sp.load_npz(f"{base}.var_called.npz").tocsc()
        m = (rl == 1) & (al == 1)
        cols = np.where(m)[0]
        # booleanise exactly as _sv_burden_build.py does: without it `~Cm` below is a bitwise
        # NOT on integers, not a mask, and the wrong genotypes get marked missing
        P = (vp[:, cols] > 0).tocsc(); Q = (vc[:, cols] > 0).tocsc()
        c = np.asarray(P.sum(0)).ravel(); k = np.asarray(Q.sum(0)).ravel()
        kk = np.where((k >= 0.9 * N) & (c >= 5) & (c <= N - 5))[0]
        kk = kk[rng.random(kk.size) < FRAC]
        G = P[:, kk].toarray().astype(np.float32)
        Cm = Q[:, kk].toarray()
        G[~Cm] = np.nan
        blocks.append(G)
        ids += [(C, int(p)) for p in pos[cols][kk]]
        print(f"  {C}: {len(kk):,} kinship SNPs", flush=True)
    return FO.founder.to_numpy(), np.hstack(blocks), ids


def grm_ours(G):
    mu = np.nanmean(G, 0); sd = np.nanstd(G, 0)
    Z = (np.where(np.isnan(G), mu, G) - mu) / sd
    return Z @ Z.T / Z.shape[1]


def write_vcf(founders, G, ids, path):
    """Inbred accessions: carrier -> 1/1, reference -> 0/0, no call -> ./."""
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in lib.CHROMS:
            f.write(f"##contig=<ID={c}>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(founders) + "\n")
        for j, (c, p) in enumerate(ids):
            g = G[:, j]
            gt = np.where(np.isnan(g), "./.", np.where(g > 0.5, "1/1", "0/0"))
            f.write(f"{c}\t{p}\t{c}_{p}\tA\tC\t.\tPASS\t.\tGT\t" + "\t".join(gt) + "\n")


def main():
    founders, G, ids = collect_snps()
    print(f"kinship genotype matrix: {G.shape[0]} founders x {G.shape[1]:,} SNPs")
    K_ours = grm_ours(G)
    np.save(f"{CHK}/grm_ours.npy", K_ours)
    K_prod = np.load(f"{OUT}/founder_kinship.npy")
    d = np.abs(K_ours - K_prod).max()
    print(f"reproduces the production matrix: max|diff| = {d:.2e}")
    assert d < 1e-5, "reconstruction does not match the production GRM -- fix before comparing"

    vcf = f"{CHK}/kinship_snps.vcf"
    write_vcf(founders, G, ids, vcf)
    cmd = [PLINK2, "--vcf", vcf, "--make-rel", "square", "--out", f"{CHK}/plink2",
           "--allow-extra-chr", "--double-id"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    print("\n".join(l for l in (p.stdout + p.stderr).splitlines()
                    if any(w in l for w in ("variants", "samples", "Error", "error", "written"))))
    K_pl = np.loadtxt(f"{CHK}/plink2.rel")
    ok = ~np.isnan(K_pl).all(1)                      # plink2 leaves all-missing samples as NaN
    if (~ok).any():
        print(f"plink2 could not place {int((~ok).sum())} sample(s); excluded from the comparison")
        K_pl = K_pl[np.ix_(ok, ok)]; K_ours = K_ours[np.ix_(ok, ok)]

    iu = np.triu_indices_from(K_ours, 1)
    r_off = np.corrcoef(K_ours[iu], K_pl[iu])[0, 1]
    r_diag = np.corrcoef(np.diag(K_ours), np.diag(K_pl))[0, 1]
    scale = np.polyfit(K_pl[iu], K_ours[iu], 1)[0]
    print(f"\noff-diagonal correlation ours vs plink2 : {r_off:.6f}")
    print(f"diagonal correlation                    : {r_diag:.6f}")
    print(f"slope (ours = slope x plink2)           : {scale:.4f}")
    print(f"ours   diag mean {np.diag(K_ours).mean():.4f}  off mean {K_ours[iu].mean():+.5f}")
    print(f"plink2 diag mean {np.diag(K_pl).mean():.4f}  off mean {K_pl[iu].mean():+.5f}")
    pd.DataFrame({"pair_i": iu[0], "pair_j": iu[1], "ours": K_ours[iu], "plink2": K_pl[iu]}) \
        .to_csv(f"{CHK}/comparison.csv", index=False)
    print(f"[wrote] {CHK}/comparison.csv")


if __name__ == "__main__":
    main()
