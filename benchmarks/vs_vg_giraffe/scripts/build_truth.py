#!/usr/bin/env python3
"""Independent (var_pa-free) truth AF for a g0 pool at the shared SNP sites.

Reads founder GTs from the ORIGINAL panel VCF (un-imputed, so '.' = not called),
restricts to the shared-panel SNP positions, and computes the MAR pool AF:

    truth_af[r] = (Σ_f w_f · alt_f) / (Σ_f w_f · called_f)

w_f from pool_weights.tsv (g0: each individual IS one founder, so founder weight
= realized pool fraction). MAR = denominator counts only CALLED founders, matching
how kMate global-mode projects AF and the sims' compute_recomb_truth convention —
but derived here straight from the founder VCF, NOT from kMate's var_pa matrix, so
it is a genuinely independent yardstick for every tool.

g0 ONLY (no recombination → no mosaics). Recombinant pools (g1/g3) need ancestry.tsv.
"""
import argparse, subprocess, sys
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--panel-vcf", required=True, help="original (un-imputed) founder VCF")
ap.add_argument("--sites-vcf", required=True, help="shared SNP panel VCF (site set)")
ap.add_argument("--weights", required=True, help="pool_weights.tsv (founder, count, weight)")
ap.add_argument("--out", required=True)
a = ap.parse_args()

# founder -> weight, aligned to panel VCF sample column order
samples = subprocess.run(["bcftools", "query", "-l", a.panel_vcf],
                         capture_output=True, text=True, check=True).stdout.split()
pw = pd.read_csv(a.weights, sep="\t", dtype={"founder": str})
wmap = dict(zip(pw["founder"], pw["weight"]))
w = np.array([wmap.get(s, 0.0) for s in samples], dtype=np.float64)
print(f"[truth] {len(samples)} panel samples; {int((w>0).sum())} with weight; Σw={w.sum():.4f}",
      file=sys.stderr)

# shared SNP positions (one biallelic SNP per pos)
sp = subprocess.run(f"bcftools view -H {a.sites_vcf} | cut -f2",
                    shell=True, capture_output=True, text=True, check=True)
shared = set(int(x) for x in sp.stdout.split())
print(f"[truth] {len(shared)} shared SNP positions", file=sys.stderr)

# stream original biallelic SNPs; compute MAR truth at shared positions
proc = subprocess.Popen(
    f"bcftools view -H -v snps -m2 -M2 {a.panel_vcf}",
    shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
rows = []
for line in proc.stdout:
    f = line.rstrip("\n").split("\t")
    pos = int(f[1])
    if pos not in shared:
        continue
    gts = f[9:]
    alt = np.zeros(len(gts)); called = np.zeros(len(gts))
    for i, g in enumerate(gts):
        c = g[0]
        if c == "0":
            called[i] = 1.0
        elif c == "1":
            called[i] = 1.0; alt[i] = 1.0
        # '.' or other -> not called
    den = float((w * called).sum())       # MAR denominator (called founders only)
    num = float((w * alt).sum())          # alt-carrying weight
    wall = float(w.sum())                 # physical denominator (ALL founders)
    rows.append((f[0], pos, f[3], f[4],
                 num / den if den > 0 else np.nan,   # truth_af   (MAR)
                 num / wall if wall > 0 else np.nan,  # truth_af_phys (missing->REF)
                 int(called.sum()), den))
proc.wait()

# truth_af      = MAR: AF among CALLED founders (matches kMate's var_pa convention)
# truth_af_phys = physical pool AF: missing GT counts as REF, because the founder
#                 FASTAs (bcftools consensus -H 1) put REF at missing sites. This is
#                 what is actually in the simulated reads -> the fair ground truth for
#                 every read-based estimator (hapFIRE, vg). Use as the primary yardstick.
df = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt",
                                 "truth_af", "truth_af_phys", "n_called", "w_called"])
df.to_csv(a.out, sep="\t", index=False)
print(f"[truth] wrote {len(df)} rows -> {a.out}", file=sys.stderr)
