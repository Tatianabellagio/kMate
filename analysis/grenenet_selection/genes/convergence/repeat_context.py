"""Repeat signatures for REF-vs-ALT motif turnover (used by tfbs_pool.py --resummarize).

tfbs_turnover.repeat_tract recognises only perfect (CT)n / (GA)n-style runs in the matched
site. Two failure modes pass it, and both inflate the 'sites lost / gained' count:

1. Redundant matching. AT4G11800's 800 bp deletion carries a degenerate GTG/GATG-rich tract
   (GGTGGTGATGATGGTGGTG...). ERF and GATA matrices match it at nearly every offset: 307
   sites lost, but only 27 distinct motifs from 3 families -- 11.4 sites per motif, against a
   pool median of 1.0 and p95 of 3.6. Counting DISTINCT motifs removes the inflation;
   sites-per-motif flags it.
2. Microsatellite length variation. A 1-17 bp indel inside a tandem repeat shifts repeat
   units, so FIMO sees sites appear and disappear that are really the same repeat. Caught by
   tandem period identity of the allele +/-25 bp: the best fraction of positions i with
   s[i] == s[i+p] (p = 1..12) in any 50 bp window.

Calibration (147 variants with turnover, 2026-09-19). Neither DUST triplet complexity nor
tandem identity over the whole allele separates AT4G11800 from GPX6 -- GPX6's 1.2 kb deletion
has an AT-rich stretch just as low-complexity (DUST >= 3 over 58% of the allele, vs 0% for
AT4G11800). Only redundancy does. Tandem identity >= 0.5 marks 13 variants, all short indels
(<= 27 bp) in microsatellites.
"""
import numpy as np

SPM_MAX = 4.0          # sites per distinct motif at or above which a variant is flagged
TANDEM_MIN = 0.5       # fraction of the allele (+flanks) lying in a tandem-periodic window


def tandem_identity(s: str, pmax: int = 12) -> float:
    s = s.upper()
    return max(np.mean([s[i] == s[i + p] for i in range(len(s) - p)])
               for p in range(1, min(pmax, len(s) - 2) + 1))


def tandem_fraction(s: str, w: int = 50, thr: float = 0.65) -> float:
    """Fraction of `s` covered by w-bp windows whose tandem identity is >= thr."""
    if len(s) < w:
        return float(tandem_identity(s) >= thr)
    cov = np.zeros(len(s), bool)
    for i in range(0, len(s) - w + 1, 5):
        if tandem_identity(s[i:i + w]) >= thr:
            cov[i:i + w] = True
    return float(cov.mean())
