#!/usr/bin/env python
"""Is the per-garden GWAS hit list driven by low-frequency markers, and is that stratum
calibrated?

Motivation (2026-08-25): the GEMMA per-garden scan returns Bonferroni hits that are
overwhelmingly low-MAC -- 90.2% of significant SNP markers sit in MAC 5-11, a stratum that is
only 36.6% of the panel, and the hit rate falls ~23-fold from rare to common. That gradient runs
BACKWARDS from statistical power: at fixed effect size a common variant is easier to detect, not
harder. So either rare variants genuinely carry much larger effects here, or the low-MAC tail is
miscalibrated. A genome-wide lambda near 1.0 cannot distinguish these -- lambda is a MEDIAN
statistic, dominated by the null bulk, and says nothing about the extreme tail where hits live.

This script stratifies by MAC and, within each stratum, reports:
  * lambda (median chi2 / 0.4549) -- bulk calibration within the stratum
  * observed/expected counts at p<1e-4 and p<1e-6 -- TAIL calibration, which is what matters
Under a well-calibrated null, obs/exp ~ 1 in every stratum. A stratum where lambda ~ 1 but the
tail obs/exp is large is inflated exactly where the hit list is drawn from.

Uses the stored Wald z directly (chi2 = z^2) rather than inverting p -- scipy's chi2.isf on
~52M p-values is the difference between seconds and >10 minutes.

Reads results/gemma_gwas/persite_gwas_{class}.npz + the arch3 panel (for MAC, cached).
Writes results/gemma_gwas/mac_calibration.csv. Env: kmate. Compute node.
"""
from __future__ import annotations
import os, sys
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from class_split_gwas import (_panel_order, _load_chrom_212, _classes, CHROMS,
                              MIN_MAC_TEST, CALL_MIN)

MIN_MAC = int(os.environ.get("KMATE_MIN_MAC", MIN_MAC_TEST))
TAG = "" if MIN_MAC == MIN_MAC_TEST else f"_mac{MIN_MAC}"
OUT = f"{lib.GEA}/r3_persite_gwas/results/gemma_gwas{TAG}"
CLASSES = ("snp", "nonsnp", "sv")
_ALL_BINS = [(5, 7), (8, 11), (12, 22), (23, 45), (46, 115)]
# drop strata the current MAC floor excludes, so the table has no empty rows
BINS = [(lo, hi) for lo, hi in _ALL_BINS if hi >= MIN_MAC]
BINS = [(max(lo, MIN_MAC), hi) for lo, hi in BINS]
CHI2_MED = stats.chi2.ppf(0.5, 1)


def mac_arrays():
    """Per-class MAC in the SAME order the scan/assemble concatenated. Cached."""
    cache = f"{OUT}/mac_by_class.npz"
    if os.path.exists(cache):
        z = np.load(cache)
        return {c: z[c] for c in CLASSES}
    order, _, _, _ = _panel_order(); N = len(order)
    acc = {c: [] for c in CLASSES}
    for cl in CHROMS:
        pos, rl, al, vp, vc, n_alt, n_cal = _load_chrom_212(cl, order)
        mac = np.minimum(n_alt, N - n_alt); vcl = _classes(rl, al)
        keep = (mac >= MIN_MAC) & (n_cal >= CALL_MIN * N)   # MUST match the scan's floor
        for name, m in (("snp", vcl == 0), ("nonsnp", vcl != 0), ("sv", vcl == 2)):
            acc[name].append(mac[keep & m])
        print(f"  MAC {cl} done", flush=True)
    out = {c: np.concatenate(v) for c, v in acc.items()}
    np.savez(cache, **out)
    return out


def main():
    mac_by = mac_arrays()
    rows = []
    print("\nTail calibration within MAC strata, pooled over all 30 gardens.")
    print("obs/exp at p<1e-4 and p<1e-6 is the TAIL check; lambda is the bulk check.\n")
    print(f"{'class':7s} {'stratum':11s} {'markers':>9s} {'lambda':>7s} "
          f"{'obs/exp p<1e-4':>16s} {'obs/exp p<1e-6':>16s}")
    for cls in CLASSES:
        z = np.load(f"{OUT}/persite_gwas_{cls}.npz", allow_pickle=True)
        Z, mac = z["Z"], mac_by[cls]
        assert len(mac) == Z.shape[0], f"{cls}: MAC/Z misalignment {len(mac)} vs {Z.shape[0]}"
        # thresholds on chi2 = z^2, avoiding a 52M-element chi2.isf
        t4 = stats.chi2.isf(1e-4, 1); t6 = stats.chi2.isf(1e-6, 1)
        for lo, hi in BINS:
            m = (mac >= lo) & (mac <= hi)
            if m.sum() < 100:
                continue
            chi = (Z[m] ** 2).ravel()
            chi = chi[np.isfinite(chi)]
            n = len(chi)
            lam = float(np.median(chi) / CHI2_MED)
            o4, e4 = int((chi > t4).sum()), n * 1e-4
            o6, e6 = int((chi > t6).sum()), n * 1e-6
            r4, r6 = o4 / e4, o6 / max(e6, 1e-12)
            lab = f"MAC {lo}-{hi}" if hi < 115 else f"MAC>={lo}"
            print(f"{cls:7s} {lab:11s} {int(m.sum()):9,d} {lam:7.3f} "
                  f"{o4:7d}/{e4:8.1f}={r4:5.2f}x {o6:5d}/{e6:7.2f}={r6:6.2f}x")
            rows.append(dict(cls=cls, mac_lo=lo, mac_hi=hi, n_markers=int(m.sum()),
                             n_tests=n, lambda_gc=lam, obs_1e4=o4, exp_1e4=e4, ratio_1e4=r4,
                             obs_1e6=o6, exp_1e6=e6, ratio_1e6=r6))
        print()
    import csv
    with open(f"{OUT}/mac_calibration.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"wrote {OUT}/mac_calibration.csv")


if __name__ == "__main__":
    main()
