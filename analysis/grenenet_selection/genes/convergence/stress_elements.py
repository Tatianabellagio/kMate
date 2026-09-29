#!/usr/bin/env python
"""Does a promoter deletion remove stress-response cis-elements -- more than a promoter
window normally carries?

1001T is rosettes in one condition, so a variant that changes a gene's INDUCIBILITY (cold,
drought, oxidative stress) rather than its baseline level is invisible there. The sequence
can still be asked: are canonical stress elements over-represented in the deleted DNA?

Elements (IUPAC, counted on both strands, non-overlapping per strand):
  CRT/DRE   RCCGAC     C-repeat / dehydration-responsive element, bound by CBF/DREB1
                       (cold) and DREB2 (drought)
  LTRE      CCGAAA     low-temperature-responsive element
  ABRE      ACGTGK     ABA-responsive element (drought, osmotic)
  W-box     TTGACY     WRKY binding (biotic / abiotic stress)
  as-1      TGACG      TGA-factor element (oxidative stress, salicylic acid)
  HSE       GAANNTTC   heat-shock element

Background: the SAME window relative to the TSS (same upstream offsets, same length) for
every TAIR10 protein-coding gene, so position and length are matched exactly. Reported per
element: count in the deleted DNA, background mean, and the percentile.

Also scans each candidate's deletion; default is the promoter deletions of the shortlist.
Writes results/stress_elements.csv.
env: kmate. Compute node, a few minutes.
"""
import os, sys, re
import numpy as np, pandas as pd, pysam
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", ".."))); sys.path.insert(0, HERE)
import lib                                                       # noqa: E402
import tfbs_turnover as TF                                       # noqa: E402

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "K": "[GT]",
         "M": "[AC]", "S": "[CG]", "W": "[AT]", "N": "[ACGT]"}
ELEM = {"CRT/DRE": "RCCGAC", "LTRE": "CCGAAA", "ABRE": "ACGTGK", "W-box": "TTGACY",
        "as-1": "TGACG", "HSE": "GAANNTTC"}
COMP = str.maketrans("ACGTRYKMSWN", "TGCAYRMKSWN")
PAT = {k: (re.compile("".join(IUPAC[c] for c in v)),
           re.compile("".join(IUPAC[c] for c in v.translate(COMP)[::-1]))) for k, v in ELEM.items()}
CAND = {"GPX6": ("AT4G11600", "Chr4", 7011705, 1225), "CRK14": ("AT4G23220", "Chr4", 12157244, 64),
        "PGIP1": ("AT5G06860", "Chr5", 2131980, 113), "CYP28": ("AT5G35100", "Chr5", 13362159, 864),
        "AT3G43230": ("AT3G43230", "Chr3", 15205485, 1580), "PKP3": ("AT1G32440", "Chr1", 11710599, 1275)}


def count(seq):
    seq = seq.upper()
    out = {}
    for k, (f, r) in PAT.items():
        n = len(f.findall(seq))
        if r.pattern != f.pattern:          # palindromes (ABRE core, as-1) not double-counted
            n += len(r.findall(seq))
        out[k] = n
    return out


def window(fa, ch, strand, tss, up_near, up_far):
    """REF sequence from up_far to up_near bp upstream of the TSS (gene orientation)."""
    if strand == "+":
        lo, hi = tss - up_far, tss - up_near
    else:
        lo, hi = tss + up_near, tss + up_far
    if lo < 1:
        return None
    return fa.fetch(ch, lo - 1, hi)


def main():
    fa = pysam.FastaFile(TF.REF_FA)
    G = lib.load_genes()
    G = G[G.chrom.isin([f"Chr{i}" for i in range(1, 6)])]
    G["tss"] = np.where(G.strand == "+", G.start, G.end)
    rows = []
    for sym, (gene, ch, pos, rl) in CAND.items():
        g = G[G.gene == gene].iloc[0]
        a, b = pos, pos + rl - 1                                  # deleted REF footprint
        if g.strand == "-":
            up_near, up_far = a - g.tss, b - g.tss
        else:
            up_near, up_far = g.tss - b, g.tss - a
        if up_near < 0:
            print(f"  {sym}: footprint reaches past the TSS -- skipped"); continue
        obs = count(window(fa, ch, g.strand, g.tss, up_near, up_far))
        bg = []
        for r in G.itertuples():
            s = window(fa, r.chrom, r.strand, r.tss, up_near, up_far)
            if s and s.upper().count("N") < 0.1 * len(s):
                bg.append(count(s))
        B = pd.DataFrame(bg)
        for k in ELEM:
            rows.append(dict(symbol=sym, element=k, window=f"-{up_far}..-{up_near}",
                             length=up_far - up_near + 1, observed=obs[k],
                             bg_mean=round(B[k].mean(), 2),
                             pctile=round(100 * (B[k] < obs[k]).mean() + 50 * (B[k] == obs[k]).mean(), 1),
                             frac_bg_ge=round((B[k] >= obs[k]).mean(), 4), n_bg=len(B)))
        print(f"  {sym:<10} window -{up_far}..-{up_near} ({up_far-up_near+1} bp) vs {len(B):,} genes")
    R = pd.DataFrame(rows)
    R.to_csv(f"{HERE}/results/stress_elements.csv", index=False)
    print(R.pivot(index="symbol", columns="element", values="observed").to_string())


if __name__ == "__main__":
    main()
