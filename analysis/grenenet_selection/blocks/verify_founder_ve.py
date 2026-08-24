#!/usr/bin/env python
"""(1) sanity-check pc1_var_explained on known inputs; (2) compute FOUNDER PC1-VE
per block and compare to the evolved-pool PC1-VE — founder-high + evolved-low = the
recombination signature. Chr1, CLQcut 0.9. Run from repo root, kmate env,
PYTHONPATH=analysis/grenenet_selection."""
import numpy as np, pandas as pd
from eval_block_coherence import pc1_var_explained
from recompute_blocks import build_common_matrix

# ---- (1) sanity tests ----
rng = np.random.RandomState(0)
F = 231
perfect = np.repeat(rng.randn(F,1), 5, axis=1)              # 5 identical cols -> VE=1
indep   = rng.randn(F, 5)                                   # independent -> VE~1/5=0.2
two_hap  = np.hstack([np.repeat(rng.randn(F,1),3,axis=1),   # 2 haplotypes (3+3 cols) -> VE~0.5
                      np.repeat(rng.randn(F,1),3,axis=1)])
print("SANITY: perfect-LD VE=%.3f (expect ~1) | independent VE=%.3f (expect ~0.2) | 2-haplotype VE=%.3f (expect ~0.5)"
      % (pc1_var_explained(perfect), pc1_var_explained(indep), pc1_var_explained(two_hap)))

# ---- (2) founder PC1-VE for CLQcut 0.9 blocks ----
print("\nbuilding founder common matrix (Chr1) ...")
std, raw, positions = build_common_matrix("panel/arch3/chr1/var_pa_231_arch3_chr1", maf=0.05, min_called_frac=0.5)
positions = np.asarray(positions)
print(f"  founders x common-variants: {std.shape}")

bl = pd.read_csv("analysis/grenenet_selection/blocks/results/blocks_recompute/chr1_clq0.9_blocks_clq0.9.tsv", sep="\t")
fves = []
for s, e in zip(bl.start_pos, bl.end_pos):
    lo, hi = np.searchsorted(positions, s), np.searchsorted(positions, e, side="right")
    fves.append(pc1_var_explained(std[:, lo:hi]) if hi-lo >= 2 else np.nan)
bl["founder_ve"] = fves

# merge evolved-VE (from the notebook data) by start_pos
ev = pd.read_csv("analysis/grenenet_selection/blocks/results/blocks_recompute/chr1_pc1ve_by_clqcut.csv")
ev = ev[ev.clqcut == 0.9][["start_pos","pc1_ve"]].rename(columns={"pc1_ve":"evolved_ve"})
m = bl.merge(ev, on="start_pos", how="left").dropna(subset=["founder_ve","evolved_ve"])

print(f"\nCLQcut 0.9, {len(m)} blocks with both VEs:")
print(f"  median FOUNDER  PC1-VE: {m.founder_ve.median():.3f}  (by-construction LD coherence)")
print(f"  median EVOLVED  PC1-VE: {m.evolved_ve.median():.3f}  (gen9 pools)")
print(f"  median drop (founder - evolved): {(m.founder_ve - m.evolved_ve).median():.3f}")
# recombination signature: tight in founders, decorrelated in evolved
recomb = m[(m.founder_ve >= 0.9) & (m.evolved_ve < 0.5)]
print(f"  blocks tight-in-founders(>=0.9) BUT low-in-evolved(<0.5) = recombination candidates: "
      f"{len(recomb)} ({100*len(recomb)/len(m):.1f}%)")
# is the low evolved tail explained by founders being low too (=block ill-defined)?
lowev = m[m.evolved_ve < 0.5]
print(f"  among low-evolved(<0.5) blocks: median founder_ve = {lowev.founder_ve.median():.3f} "
      f"(if high -> NOT a block-definition problem -> recombination/noise)")
m.to_csv("analysis/grenenet_selection/blocks/results/blocks_recompute/chr1_clq0.9_founder_vs_evolved_ve.csv", index=False)
print("wrote chr1_clq0.9_founder_vs_evolved_ve.csv")
