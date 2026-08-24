#!/usr/bin/env python
"""WZA core, parameterized for the wza_investigation.

Faithful re-implementation of Booker et al. 2024 (general_WZA_script.py) but with
the SNP-number-correction degree, the rolling-window params, and the per-window
SNP cap all exposed so we can sweep them. The raw weighted-Z is byte-identical to
the canonical script; only the correction stage is parameterized.

raw block stat:   Z = Σ(pq·z) / sqrt(Σ pq²),   pq = MAF(1-MAF),  z = Φ⁻¹(1-p)
optional cap:     if n_snps > cap, average Z over `resamples` random cap-subsets
correction:       sort by SNPs; rolling(roller, min=minEntries) mean/var of Z;
                  polyfit(mean_SNP, {mean,sd}, deg); Z_pVal = 1-Φ(Z; mean̂, sd̂)
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import scipy.stats
from scipy.stats import norm

RNG = np.random.default_rng(0)  # fixed seed -> reproducible cap resampling


def _block_Z(sub: pd.DataFrame) -> float:
    """Raw weighted-Z for one window (sub has columns z_score, pbar_qbar)."""
    num = (sub["pbar_qbar"] * sub["z_score"]).sum()
    den = np.sqrt((sub["pbar_qbar"] ** 2).sum())
    return np.nan if den == 0 else num / den


def raw_wza(df: pd.DataFrame, pcol="pval", maf="MAF", window="block",
            cap=None, resamples=100, min_snps=2) -> pd.DataFrame:
    """One raw weighted-Z per window. `cap` (int) downsamples big windows with
    `resamples` averaging, exactly as Booker --sample_snps. cap=None -> no cap."""
    d = df.copy()
    p = d[pcol].clip(lower=1e-15).replace(1, 1 - 1e-3).astype(float)
    d["z_score"] = scipy.stats.norm.ppf(1 - p.to_numpy())
    d["pbar_qbar"] = d[maf] * (1 - d[maf])
    rows = []
    for w, sub in d.groupby(window):
        n = len(sub)
        if n < min_snps:
            continue
        if cap is None or n <= cap:
            Z = _block_Z(sub); used = n
        else:
            Z = np.mean([_block_Z(sub.iloc[RNG.choice(n, cap, replace=False)])
                         for _ in range(resamples)])
            used = cap
        rows.append((w, n, used, Z, sub[maf].mean()))
    return pd.DataFrame(rows, columns=["block", "SNPs_raw", "SNPs", "Z", "MAF"])


def fit_correction(wza: pd.DataFrame, deg=2, roller=50, minEntries=40):
    """Return (mean_poly, sd_poly, diag) fit on rolling bins of Z vs SNP count.
    diag carries the rolling support points so we can plot the fit + the tail."""
    s = wza[~wza.Z.isnull()].sort_values("SNPs")
    var = s.Z.rolling(roller, min_periods=minEntries).var()
    mask = ~var.isnull()
    sd = np.sqrt(var[mask])
    mean = s.Z.rolling(roller, min_periods=minEntries).mean()[mask]
    xsnp = s.SNPs.rolling(roller, min_periods=minEntries).mean()[mask]
    sd_poly = np.poly1d(np.polyfit(xsnp, sd, deg))
    mean_poly = np.poly1d(np.polyfit(xsnp, mean, deg))
    diag = dict(x=xsnp.to_numpy(), sd=sd.to_numpy(), mean=mean.to_numpy())
    return mean_poly, sd_poly, diag


def apply_correction(wza: pd.DataFrame, deg=2, roller=50, minEntries=40,
                     sd_floor=False):
    """Add Z_pVal. NO safety floor by default -> exposes the NaN failure mode.
    sd_floor=True clips predicted SD to smallest positive (the kMate hack)."""
    wza = wza.copy()
    mean_poly, sd_poly, diag = fit_correction(wza, deg, roller, minEntries)
    sd_pred = sd_poly(wza["SNPs"].to_numpy())
    mean_pred = mean_poly(wza["SNPs"].to_numpy())
    wza["sd_pred"] = sd_pred
    wza["mean_pred"] = mean_pred
    wza["neg_sd"] = sd_pred <= 0
    if sd_floor:
        pos = sd_pred[sd_pred > 0]
        floor = pos.min() if pos.size else np.nanstd(wza["Z"].to_numpy())
        sd_pred = np.clip(sd_pred, floor, None)
    with np.errstate(invalid="ignore"):
        wza["Z_pVal"] = 1 - norm.cdf(wza["Z"].to_numpy(), loc=mean_pred, scale=sd_pred)
    return wza, diag
