#!/usr/bin/env python3
"""Panel / filter comparison of per-block AF error (w10kb, floor=1, no smoothing).

Three series on the SAME 3 scenarios (n50_g0, n231_g1_self97, n50_g3_dom500nr):
  p231 filt2inv  : 231 mixed-tech founders, private (ac=1) DROPPED   (production)
  p80  filt2     : 80 all-long-read founders, private DROPPED
  p80  unfiltered: 80 all-long-read founders, private KEPT           ("delete the filter")

Answers: (1) does the 80-founder long-read panel resolve better? [confounded by
fewer founders]; (2) CLEAN within-p80 test — does keeping private k-mers help?
Plus a centromere-split, since that's where panel k-mer quality should matter most."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = Path("benchmarks/localonly_p231/results")
KEYS = ["chrom", "pos", "ref_len", "alt_len"]
SCEN = ["cov10_n50_g0_s42_hotspots", "cov10_n231_g1_s42_self97_hotspots",
        "cov10_n50_g3_s42_hotspots_dom500nr"]
EDGES = [1, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 100000]
# (label, color, floor_diag dir, pool suffix, variant tag, sims dir, truth filename)
SERIES = [
    ("p231 filt2inv (231 mixed)", "C0", "floor_diag", "_p231_chr1", "",
     "benchmarks/p231/sims", "recomb_truth_raw.tsv.gz"),
    ("p80 filt2 (80 long-read)", "C1", "floor_diag_p80", "_p80_chr1", "filt2",
     "benchmarks/p80/sims", "recomb_truth.tsv.gz"),
    ("p80 unfiltered (keep private)", "C2", "floor_diag_p80", "_p80_chr1", "unfiltered",
     "benchmarks/p80/sims", "recomb_truth.tsv.gz"),
]
CEN = (12_500_000, 17_500_000)


def per_block(fd, pool, variant, sims, truthfn, min_rec=20):
    tag = f"{pool}_w10kb" + (f"_{variant}" if variant else "")
    bd = Path("benchmarks/localonly_p231") / fd / f"{tag}.blockdiag.tsv"
    ef = Path("benchmarks/localonly_p231") / fd / f"{tag}.recest.tsv"
    if not bd.exists() or not ef.exists():
        return None
    diag = pd.read_csv(bd, sep="\t"); est = pd.read_csv(ef, sep="\t")
    tr = pd.read_csv(Path(sims) / pool / truthfn, sep="\t").dropna(subset=["truth_af"])
    tr["occ"] = tr.groupby(KEYS).cumcount(); est["occ"] = est.groupby(KEYS).cumcount()
    m = tr.merge(est[KEYS + ["occ", "alt_freq", "block"]], on=KEYS + ["occ"], how="inner")
    m = m.merge(diag[["block", "nnz", "start"]], on="block", how="left")
    m = m[np.isfinite(m.alt_freq.values)].copy()
    m["e2"] = (m.alt_freq - m.truth_af) ** 2
    g = pd.DataFrame({"rmse": np.sqrt(m.groupby("block").e2.mean()),
                      "n": m.groupby("block").alt_freq.size(),
                      "nnz": m.groupby("block").nnz.first(),
                      "start": m.groupby("block").start.first()})
    g = g[g.n >= min_rec].copy()
    g["cen"] = (g.start >= CEN[0]) & (g.start <= CEN[1])
    return g


# gather all series
data = {}
for lab, col, fd, suf, var, sims, tfn in SERIES:
    parts = [per_block(fd, scn + suf, var, sims, tfn) for scn in SCEN]
    parts = [p for p in parts if p is not None]
    if parts:
        data[lab] = (col, pd.concat(parts, ignore_index=True))
        print(f"loaded {lab}: {len(data[lab][1])} blocks")
    else:
        print(f"MISSING {lab} (job not done?)")

if not data:
    raise SystemExit("no inputs yet")


def med_curve(df):
    b = df.copy(); b["bin"] = pd.cut(b.nnz, EDGES, right=False)
    ctr = b.groupby("bin", observed=True).nnz.median()
    med = b.groupby("bin", observed=True).rmse.median()
    return ctr.values, med.values


fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
for title, sub in [("ALL blocks", lambda d: d), ("CENTROMERE blocks only", lambda d: d[d.cen])]:
    ax = axes[0] if title == "ALL blocks" else axes[1]
    for lab, (col, df) in data.items():
        d = sub(df)
        if len(d) < 30:
            continue
        x, y = med_curve(d)
        ax.plot(x, y, "-o", color=col, label=lab, ms=4)
    ax.axhline(0.048, color="gray", ls=":", lw=1)
    ax.set_xscale("log"); ax.set_xlabel("nonzero k-mers in block (nnz)")
    ax.set_title(title, fontsize=10); ax.set_ylim(0, 0.30)
axes[0].set_ylabel("median per-block AF RMSE"); axes[0].legend(fontsize=8)
fig.suptitle("Panel & private-k-mer-filter effect on per-block AF error (w10kb)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(RES / "panel_compare_floor.png", dpi=140)
print(f"-> {RES/'panel_compare_floor.png'}")

# numeric summary
print("\n=== centromere penalty at matched supply (nnz>=512) ===")
print(f"{'series':>32} {'arm_med':>8} {'cen_med':>8} {'penalty':>8} {'cen_med_nnz':>11}")
for lab, (col, df) in data.items():
    hi = df[df.nnz >= 512]
    arm = hi[~hi.cen].rmse.median(); cen = hi[hi.cen].rmse.median()
    print(f"{lab:>32} {arm:>8.3f} {cen:>8.3f} {cen-arm:>+8.3f} {hi[hi.cen].nnz.median():>11.0f}")

# how many k-mers does keeping private add, esp in the centromere?
if "p80 filt2 (80 long-read)" in data and "p80 unfiltered (keep private)" in data:
    f2 = data["p80 filt2 (80 long-read)"][1]; uf = data["p80 unfiltered (keep private)"][1]
    print("\n=== median block nnz: filt2 vs unfiltered (p80) ===")
    print(f"  arm:        filt2={f2[~f2.cen].nnz.median():.0f}  unfiltered={uf[~uf.cen].nnz.median():.0f}")
    print(f"  centromere: filt2={f2[f2.cen].nnz.median():.0f}  unfiltered={uf[uf.cen].nnz.median():.0f}")
