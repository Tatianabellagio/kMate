#!/usr/bin/env python
"""Raw AF-vs-climate + kMate panel support for the Chr4 CRK block (4_2519, CRK10-14).

Three questions:
  (A) Did the block's allele frequency rise in warm or cold sites? (raw gen9 AF vs bio1,
      site-level; + Delta p = gen9 - founding p0, the selection signature).
  (B) How much PANEL SUPPORT did this allele have? (kMate n_called = # of the 231 founder
      panel haplotypes informing each call; from the raw per-chrom TSVs).
  (C) The kMate ERROR (se) on these calls.

CRK10-14 = AT4G23180..AT4G23220 = Chr4:12,138,137-12,157,091, inside LD block 4_2519.

Outputs: results/.../phase1_replication/crk_block_af_support.png + printed stats.
  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python (basic env for plotting)
"""
from __future__ import annotations
import glob, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import kendalltau
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))  # followups/ -> grenenet_gea
import lib

CM = f"{lib.GEA}/phase1_replication/class_matrices"
STORE = lib.AF_STORE
BASE = "/global/scratch/users/tbellg/kmate/results/grenenet_kmate_arch3"
CHROM, LO, HI = "Chr4", 12138137, 12157091     # CRK10-14 span
BLOCK = "4_2519"
OUT = f"{lib.GEA}/phase1_replication/followups/crk_block_af_support.png"


def main():
    # ---------- (A) raw AF vs climate, site level ----------
    recs = pd.read_csv(f"{CM}/snp_gen9.records.csv")
    pools = pd.read_csv(f"{CM}/gen9.pools.csv")
    af = np.load(f"{CM}/snp_gen9_af.npy")                       # [pools x records]
    sel = ((recs.chrom == CHROM) & (recs.pos >= LO) & (recs.pos <= HI)).to_numpy()
    if "block" in recs:
        sel &= (recs.block.astype(str) == BLOCK).to_numpy()
    idx = np.where(sel)[0]
    print(f"(A) CRK10-14 SNPs in block {BLOCK}: {len(idx)}")
    A = af[:, idx]                                              # [pools x crkSNP]
    bio1 = pools["bio1"].to_numpy(float)
    site = pools["site"].to_numpy()
    fw = pools["total_flowers"].to_numpy(float)

    # founding p0 for these SNPs (map by chrom,pos to the af_store index)
    ix = np.load(f"{STORE}/index_snp.npz", allow_pickle=True)
    p0all = np.load(f"{STORE}/p0_snp.npy")
    key = {(c, int(p)): i for i, (c, p) in enumerate(zip(ix["chrom"], ix["pos"]))}
    p0 = np.array([p0all[key.get((CHROM, int(p)), -1)] for p in recs.pos.to_numpy()[idx]])

    # per-SNP spatial Kendall tau(ALT-AF, bio1) across pools -> direction
    taus = []
    for j in range(A.shape[1]):
        y = A[:, j]; ok = np.isfinite(y)
        if ok.sum() >= 3 and np.ptp(y[ok]) > 0:
            taus.append(kendalltau(bio1[ok], y[ok])[0])
    taus = np.array(taus)
    print(f"    per-SNP tau(ALT-AF, bio1): mean {taus.mean():+.3f} | frac<0 {(taus<0).mean():.2f} "
          f"-> ALT allele {'DOWN in warm / UP in cold' if taus.mean()<0 else 'UP in warm'}")

    # orient every SNP to its cold-favoured allele (flip ALT where tau>0) so the block
    # haplotype frequency is a single coherent cline (visualization only; tau test above is raw)
    orient = np.where(taus[None, :] > 0, 1.0 - A, A)
    # site-level flower-weighted mean of the oriented block frequency, gen9 + delta-p
    df = pd.DataFrame({"site": site, "bio1": bio1, "fw": fw,
                       "hap": np.nanmean(orient, axis=1)})
    p0_oriented = np.where(taus > 0, 1.0 - p0, p0)
    df["dp"] = np.nanmean(orient - p0_oriented[None, :], axis=1)
    g = (df.groupby("site")
           .apply(lambda d: pd.Series({
               "bio1": d.bio1.iloc[0],
               "hap": np.average(d.hap, weights=d.fw),
               "dp": np.average(d.dp, weights=d.fw)}))
           .reset_index())
    tg = kendalltau(g.bio1, g.hap); td = kendalltau(g.bio1, g.dp)
    print(f"    SITE-level (n={len(g)}): hap-freq vs bio1 tau={tg[0]:+.3f} p={tg[1]:.2g} | "
          f"deltap vs bio1 tau={td[0]:+.3f} p={td[1]:.2g}")

    # ---------- (B/C) panel support n_called + error se from raw Chr4 TSVs ----------
    tsvs = sorted(glob.glob(f"{BASE}/*_Chr4.tsv"))[:20]
    print(f"\n(B/C) reading {len(tsvs)} sample Chr4 TSVs for n_called + se …")
    nc_crk, se_crk, nc_bg = [], [], []
    for t in tsvs:
        d = pd.read_csv(t, sep="\t", usecols=["pos", "n_called", "se"])
        c = d[(d.pos >= LO) & (d.pos <= HI)]
        nc_crk.append(c.n_called.to_numpy()); se_crk.append(c.se.to_numpy())
        nc_bg.append(d.n_called.sample(min(5000, len(d)), random_state=1).to_numpy())
    nc_crk = np.concatenate(nc_crk); se_crk = np.concatenate(se_crk)
    nc_bg = np.concatenate(nc_bg)
    print(f"    panel support n_called (CRK10-14): median {np.median(nc_crk):.0f} "
          f"mean {nc_crk.mean():.1f} [{nc_crk.min():.0f}-{nc_crk.max():.0f}]  (panel max=231)")
    print(f"    genome-wide baseline n_called:     median {np.median(nc_bg):.0f} mean {nc_bg.mean():.1f}")
    print(f"    kMate se (CRK10-14):  median {np.nanmedian(se_crk):.4f} "
          f"mean {np.nanmean(se_crk):.4f} [{np.nanmin(se_crk):.3f}-{np.nanmax(se_crk):.3f}]")

    # ---------- figure ----------
    fig, ax = plt.subplots(1, 4, figsize=(20, 4.4))
    ax[0].scatter(g.bio1, g.hap, c=g.bio1, cmap="coolwarm", s=60, edgecolor="k")
    ax[0].set(xlabel="bio1 (annual mean T, °C)", ylabel="CRK-block hap freq (gen9)",
              title=f"(A) raw AF vs climate\nsite tau={tg[0]:+.2f} p={tg[1]:.1g}")
    ax[1].scatter(g.bio1, g.dp, c=g.bio1, cmap="coolwarm", s=60, edgecolor="k")
    ax[1].axhline(0, color="grey", lw=.8)
    ax[1].set(xlabel="bio1 (°C)", ylabel="Δp = gen9 − p0",
              title=f"(A) selection signature\nsite tau={td[0]:+.2f} p={td[1]:.1g}")
    ax[2].hist(nc_bg, bins=40, alpha=.5, density=True, label="genome-wide", color="grey")
    ax[2].hist(nc_crk, bins=40, alpha=.7, density=True, label="CRK10-14", color="firebrick")
    ax[2].axvline(231, color="k", ls=":", lw=.8); ax[2].set(
        xlabel="kMate n_called (panel support)", ylabel="density",
        title=f"(B) panel support\nCRK median={np.median(nc_crk):.0f}/231"); ax[2].legend(fontsize=8)
    ax[3].hist(se_crk[np.isfinite(se_crk)], bins=40, color="steelblue")
    ax[3].set(xlabel="kMate se (AF std error)", ylabel="count",
              title=f"(C) kMate error\nmedian se={np.nanmedian(se_crk):.3f}")
    fig.suptitle("Chr4 CRK block 4_2519 (CRK10-14, Chr4:12.14-12.16 Mb) — AF vs climate + kMate support",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT, dpi=130, bbox_inches="tight")
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
