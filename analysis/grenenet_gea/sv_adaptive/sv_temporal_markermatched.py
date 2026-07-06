#!/usr/bin/env python3
"""Does the POOL-TEMPORAL block-level SV enrichment survive hap-count matching?

Companion to sv_frac_markermatched.py (founder-GWAS). The pool-temporal "all-blocks" SV test
(notebook cell 3b, memory ×1.01->×1.04 after the MAC>=2 refloor) scores each block by
sel_gea = MAX |s_mean| over the block's testable HapFM haplotypes -- a max-over-markers block
aggregation, structurally identical to the founder-GWAS block score (max chi2 over hap-cluster
markers) that collapsed x1.05->x0.999 under n_markers matching (sv_joint_mechanism.py). So the
pool-temporal x1.04 could be the SAME hap-count aggregation artifact. The per-site leg
(cross_site_enrichment, median of SNP s -> NOT a max) is immune and is the clean instrument.

Tests, on the current MAC>=2 sv_landscape:
  (a) do SV-bearing blocks have more testable haplotypes than n_kept-matched SNP blocks?
      (if yes, sel_gea = max over more haps inflates for free)
  (b) reproduce the all-blocks SV-bearing sel_gea excess, n_kept-matched (~x1.04)
  (c) re-match on n_kept x n_haps -- does the excess survive the aggregation confound?
Run in `basic` from the repo root. Deterministic (seed=0). NPERM=10000.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pathlib import Path
import numpy as np, pandas as pd
import lib

OUT = Path("results/grenenet_gea/sv_adaptive")
GEADIR = Path("results/grenenet_gea/hapfreq_clq90/pipelineB_varlen")
NPERM = 10000
EDGES = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 70, 100, 150, 250, 10**9]
HMAX = 5

hg = pd.read_csv(GEADIR / "hap_gea.csv")
hg = hg[hg.covered & hg.panel_freq.between(0.05, 0.95)].copy()      # the GEA test's testable set
hg["unit"] = hg.chrom + ":" + hg.unit_start.astype(str) + "-" + hg.unit_end.astype(str)
blk = hg.groupby("unit").agg(sel_gea=("s_mean", lambda s: np.abs(s).max()),
                             n_haps=("s_mean", "size")).reset_index()

L = pd.read_csv(OUT / "sv_landscape_clq0.9.csv")
df = blk.merge(L[["block_id", "n_kept", "has_sv", "sv_frac"]],
               left_on="unit", right_on="block_id", how="inner").reset_index(drop=True)
df["kbin"] = pd.cut(df.n_kept, EDGES, right=False, labels=False)
df = df[df.kbin.notna()].reset_index(drop=True); df["kbin"] = df.kbin.astype(int)
df["hbin"] = np.minimum(df.n_haps, HMAX)
df["jbin"] = df.kbin * 10 + df.hbin
kb, jb = df.kbin.to_numpy(), df.jbin.to_numpy()
sv = np.where(df.has_sv.to_numpy() == 1)[0]
print(f"{len(df):,} blocks with a pool-temporal score | {len(sv):,} SV-bearing | "
      f"n_haps/block: SNP-only median {df.loc[df.has_sv==0,'n_haps'].median():.0f} "
      f"mean {df.loc[df.has_sv==0,'n_haps'].mean():.2f} | "
      f"SV median {df.loc[df.has_sv==1,'n_haps'].median():.0f} mean {df.loc[df.has_sv==1,'n_haps'].mean():.2f}")

# (a) hap-count excess of SV blocks at matched n_kept (the confound magnitude)
o, n, r, p = lib.matched_perm_test(df.n_haps.to_numpy(float), kb, sv, NPERM)
print(f"\n(a) n_haps/block: SV {o:.2f} vs n_kept-matched {n:.2f}  x{r:.3f}  p={p:.4f}")

# (b) all-blocks SV-bearing sel_gea excess, n_kept-matched (reproduce memory's ~x1.04)
o1, n1, r1, p1 = lib.matched_perm_test(df.sel_gea.to_numpy(), kb, sv, NPERM)
# (c) re-matched on n_kept x n_haps
o2, n2, r2, p2 = lib.matched_perm_test(df.sel_gea.to_numpy(), jb, sv, NPERM)
print(f"\n(b) all-blocks SV-bearing sel_gea (=max|s_mean|), n_kept-matched:  "
      f"SV {o1:.4f} vs {n1:.4f}  x{r1:.3f}  p={p1:.4f}")
print(f"(c) same, matched on n_kept x n_haps:                             "
      f"SV {o2:.4f} vs {n2:.4f}  x{r2:.3f}  p={p2:.4f}  "
      f"({'SURVIVES hap-count matching' if r2 > 1 and p2 < 0.05 else 'excess GONE -> hap-count aggregation artifact'})")

df.to_csv(OUT / "sv_temporal_markermatched_blocks.csv", index=False)
print(f"\n[done] -> {OUT}/sv_temporal_markermatched_blocks.csv")
