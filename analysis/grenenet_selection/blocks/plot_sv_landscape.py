#!/usr/bin/env python3
"""Visualize the per-haploblock SV landscape + the two confounds any enrichment
test must control: (b) block size, (c) genomic context. Run in `basic` env."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("analysis/grenenet_selection/sv_adaptive/results")
df = pd.read_csv(OUT / "sv_landscape_clq0.9.csv")
CEN = {"Chr1": 15.086, "Chr2": 3.607, "Chr3": 13.799, "Chr4": 3.956, "Chr5": 11.725}
CHROMS = [f"Chr{i}" for i in range(1, 6)]

fig, ax = plt.subplots(2, 2, figsize=(13, 9))

# (a) how many SVs per block — sparsity
a = ax[0, 0]
vc = df.n_sv.clip(upper=6).value_counts().sort_index()
a.bar(vc.index, vc.values, color="#4a7", edgecolor="white")
a.set_yscale("log"); a.set_xlabel("SVs in block (n_sv, capped at 6)"); a.set_ylabel("# blocks (log)")
a.set_title(f"(a) SV content is sparse — {100*df.has_sv.mean():.1f}% of blocks have any SV", loc="left", fontsize=10)
for x, y in zip(vc.index, vc.values):
    a.text(x, y, f"{y:,}", ha="center", va="bottom", fontsize=7)

# (b) THE size confound: P(has SV) rises with block variant count
b = ax[0, 1]
df["nkbin"] = pd.cut(df.n_kept, [1, 3, 5, 8, 15, 30, 60, 120, 10000],
                     labels=["2","3-4","5-7","8-14","15-29","30-59","60-119","120+"])
g = df.groupby("nkbin", observed=True).agg(p=("has_sv","mean"), n=("has_sv","size"))
b.bar(range(len(g)), 100*g.p.values, color="#c65", edgecolor="white")
b.set_xticks(range(len(g))); b.set_xticklabels(g.index, fontsize=8)
b.set_xlabel("variants in block (n_kept)"); b.set_ylabel("% blocks with >=1 SV")
b.set_title("(b) SIZE CONFOUND — bigger blocks ⇒ more likely to carry an SV", loc="left", fontsize=10)
for i,(p,n) in enumerate(zip(g.p,g.n)):
    b.text(i, 100*p, f"{100*p:.0f}%\nn={n:,}", ha="center", va="bottom", fontsize=6.5)

# (c) genomic-context confound: SV density along each chromosome, centromere marked
c = ax[1, 0]
off = 0; ticks = []
for ch in CHROMS:
    sub = df[df.chrom == ch].sort_values("mid")
    x = sub.mid.to_numpy()/1e6 + off
    c.scatter(x, sub.sv_density_kb, s=3, alpha=0.25,
              color="#345" if CHROMS.index(ch) % 2 == 0 else "#69a")
    c.axvline(CEN[ch] + off, color="crimson", ls=":", lw=1)
    ticks.append((off + x.mean() - off + sub.mid.mean()/1e6*0 + (sub.mid.max()/1e6)/2, ch))
    off += sub.mid.max()/1e6 * 1.03
c.set_ylim(0, np.percentile(df.sv_density_kb[df.sv_density_kb>0], 99))
c.set_xlabel("genome position (Mb, concatenated)"); c.set_ylabel("SV density (SV / kb)")
c.set_title("(c) CONTEXT CONFOUND — SV density spikes at (peri)centromeres (red)", loc="left", fontsize=10)

# (d) among SV-containing blocks, what fraction of their variants are SVs
d = ax[1, 1]
sf = df.loc[df.has_sv == 1, "sv_frac"]
d.hist(sf, bins=40, color="#4a7", edgecolor="white")
d.axvline(sf.median(), color="k", ls="--", lw=1, label=f"median {sf.median():.3f}")
d.set_xlabel("SV fraction of block's variants (SV-containing blocks)"); d.set_ylabel("# blocks")
d.set_title(f"(d) where SVs occur they're a small minority (n={len(sf):,} blocks)", loc="left", fontsize=10)
d.legend(fontsize=8)

fig.suptitle("Structural-variant landscape across 58,376 clq0.9 haploblocks "
             "(founder MAF≥0.05; SV = |Δlen|>50 bp)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUT / "sv_landscape.png", dpi=140)
print(f"-> {OUT/'sv_landscape.png'}")

# quantify the confounds numerically
print("\n== SIZE confound: logistic-ish trend P(has_sv) vs log10(n_kept) ==")
lr = np.corrcoef(np.log10(df.n_kept), df.has_sv)[0, 1]
print(f"  corr(log10 n_kept, has_sv) = {lr:.3f}   (must be controlled)")
print("== CONTEXT confound: SV density in pericentromere vs arm ==")
df["peri"] = df.dist_cen < 2.5e6
print(f"  pericentromere (<2.5 Mb): {100*df[df.peri].has_sv.mean():.1f}% have SV  (n={df.peri.sum():,})")
print(f"  arm:                      {100*df[~df.peri].has_sv.mean():.1f}% have SV  (n={(~df.peri).sum():,})")
