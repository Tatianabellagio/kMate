#!/usr/bin/env python
"""Shared founder-panel genotype builder + EMMAX variance-component helper.

build_genotype() constructs the founder x haploblock-cluster 0/1 matrix from the dynld-K500
haplotype membership, plus a per-marker chrom/pos/block registry. emma_reml_delta() is the
EMMAX/P3D ML estimate of delta = sigma_e^2/sigma_g^2 on the kinship eigenbasis.

These were factored out of the old founder_gwas_site.py (the 194-founder GWAS, which was
incorrect -- it dropped founders kMate estimates at ~0, but all 231 ecotypes are in the seed
mix so those are estimation artifacts, not real absence; use founder_gwas_231 / _precision).
The functions themselves are panel infrastructure and are imported by the block LD-LMM and the
231-founder GWAS scripts. Env: kmate.
"""
import numpy as np
import os
import pandas as pd
from scipy import optimize

HM = "analysis/grenenet_selection/blocks/results/blocks_mcf90/hap_membership"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]


def build_genotype(tag=None):
    """Founder x haplotype 0/1 matrix from haplotype membership + per-hap chrom/pos.

    `tag` selects the block/membership set: 'K500' (default, dynld-K500 units) or 'clq90'
    (finer r2>=0.9 LD blocks). Overridable via the MEMB_TAG env var. The GWAS genotype is
    panel haplotype membership and the trait is genome-wide ecotype h, so the k-mer-coverage
    rationale for the coarse K500 units does not apply -- clq90 gives finer localization.
    """
    tag = tag or os.environ.get("MEMB_TAG", "K500")
    Gs, chrom, start, end, block = [], [], [], [], []
    founders = None
    boff = 0                                               # global block (unit) counter
    for ch in CHROMS:
        M = np.load(f"{HM}/{ch.lower()}_hapmemb_{tag}.npz", allow_pickle=True)
        labels = M["labels"].astype(np.int64)             # (U, F) cluster of each founder
        off = M["hap_offset"].astype(np.int64)            # (U+1,)
        nh = int(off[-1]); F = labels.shape[1]; U = len(off) - 1
        if founders is None:
            founders = M["founders"].astype(str)
        gid = off[:-1][:, None] + labels                  # (U, F) global hap id within chrom
        G = np.zeros((F, nh), np.float32)
        gidT = gid.T                                       # (F, U)
        for f in range(F):
            G[f, gidT[f]] = 1.0
        Gs.append(G)
        us, ue = M["unit_start"].astype(np.int64), M["unit_end"].astype(np.int64)
        nclust = np.diff(off)
        chrom.append(np.repeat(ch, nh))
        start.append(np.repeat(us, nclust)); end.append(np.repeat(ue, nclust))
        block.append(np.repeat(np.arange(U) + boff, nclust)); boff += U
    G = np.hstack(Gs)                                      # (F, nhap_total)
    reg = pd.DataFrame({"chrom": np.concatenate(chrom),
                        "start": np.concatenate(start), "end": np.concatenate(end),
                        "block": np.concatenate(block)})
    return G, founders, reg


def emma_reml_delta(y, X, lam):
    """ML estimate of delta = sigma_e^2/sigma_g^2 on the K-eigenbasis (y,X already rotated)."""
    n = len(y)
    def negll(logdelta):
        d = np.exp(logdelta); w = 1.0 / (lam + d)
        XtWX = X.T * w @ X; XtWy = X.T * w @ y
        beta = np.linalg.solve(XtWX, XtWy)
        resid = y - X @ beta
        s2 = (w * resid**2).sum() / n
        return 0.5 * (n * np.log(2 * np.pi * s2) + np.log(lam + d).sum() + n)
    r = optimize.minimize_scalar(negll, bounds=(-8, 8), method="bounded")
    return float(np.exp(r.x))
