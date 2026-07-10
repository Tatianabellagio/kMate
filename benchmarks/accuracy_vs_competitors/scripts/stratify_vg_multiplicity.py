#!/usr/bin/env python3
"""Diagnose WHY vg-SV underperforms kMate on the g0 pool: is it driven by
multiallelic cactus snarls (many competing ALT alleles stacked at one anchor),
where short reads cannot apportion support, as opposed to vg being inaccurate?

We stratify vg's call -v AF accuracy by the allele multiplicity of the locus
(#SV ALT alleles sharing the same POS in the panel SV VCF) and by SV size, and
overlay kMate on the SAME SVs. Output: a metrics table + a 3-panel figure.

All tools join to the SAME truth on the panel SV index (svidx = row in the
SV-masked panel meta order, identical convention to score_sv.py / build_4tool).
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/global/scratch/users/tbellg/kmate")
POOL = sys.argv[1] if len(sys.argv) > 1 else "cov10_n80_g0_s42_hotspots_p80_chr1"
SVLEN = 50

PANEL_SV = ROOT / "benchmarks/accuracy_vs_competitors/work/panel_sv_p80_Chr1.vcf.gz"
VG_VCF   = ROOT / f"benchmarks/accuracy_vs_competitors/work/{POOL}_vgconstruct_sv_FIXED.vcf"
TRUTH    = ROOT / f"benchmarks/p80/sims/{POOL}/recomb_truth.tsv.gz"
META     = ROOT / "benchmarks/p80/data/var_pa_p80.meta.npz"
CALLED   = ROOT / "benchmarks/p80/data/var_pa_p80.var_called.npz"
KMATE    = ROOT / f"benchmarks/benchmark_runs/tsv/{POOL}_global.tsv"
OUTPFX   = ROOT / f"benchmarks/accuracy_vs_competitors/results/vg_multiplicity_{POOL}"


def seq_key(chrom, pos, ref, alt):
    return f"{chrom.replace('Chr','').replace('chr','')}:{pos}:{ref}:{alt}"


def metrics(truth, est):
    m = np.isfinite(truth) & np.isfinite(est)
    truth, est = truth[m], est[m]
    if len(truth) < 2:
        return dict(n=int(m.sum()), MAE=np.nan, RMSE=np.nan, R2=np.nan, r=np.nan, slope=np.nan)
    err = est - truth
    ss_res = float((err ** 2).sum())
    ss_tot = float(((truth - truth.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    r = float(np.corrcoef(truth, est)[0, 1])
    slope = float(np.polyfit(truth, est, 1)[0])
    return dict(n=int(m.sum()), MAE=float(np.abs(err).mean()),
                RMSE=float(np.sqrt((err ** 2).mean())), R2=r2, r=r, slope=slope)


# ---- 1. panel SV VCF: svidx (mask order) + POS + per-POS multiplicity --------
print("[1] reading panel SV VCF ...", file=sys.stderr)
key2idx, pos_of_idx = {}, []
op = gzip.open if str(PANEL_SV).endswith(".gz") else open
with op(PANEL_SV, "rt") as fh:
    idx = 0
    for ln in fh:
        if ln.startswith("#"):
            continue
        f = ln.split("\t")
        chrom, pos, ref, alts = f[0], int(f[1]), f[3], f[4].split(",")
        for alt in alts:
            if abs(len(alt) - len(ref)) < SVLEN:
                continue
            key2idx[seq_key(chrom, pos, ref, alt)] = idx
            pos_of_idx.append(pos)
            idx += 1
pos_of_idx = np.array(pos_of_idx)
nsv = len(pos_of_idx)
# multiplicity = how many SV alleles share this svidx's POS
posser = pd.Series(pos_of_idx)
mult = posser.map(posser.value_counts()).values
print(f"    {nsv:,} SV alleles; {len(np.unique(pos_of_idx)):,} distinct POS; "
      f"multiplicity max={mult.max()}, median={int(np.median(mult))}", file=sys.stderr)

# ---- 2. truth (physical AF) + n_called, in SV-mask order ---------------------
print("[2] reading truth + meta ...", file=sys.stderr)
tr = pd.read_csv(TRUTH, sep="\t")
m = np.load(META, allow_pickle=True)
assert len(tr) == len(m["pos"]), f"{len(tr)} != {len(m['pos'])}"
svmask = np.abs(m["alt_len"].astype(int) - m["ref_len"].astype(int)) >= SVLEN
vc = sp.load_npz(CALLED).tocsr()
F = vc.shape[0]
n_called = np.asarray(vc[:, svmask].sum(axis=0)).ravel().astype(int)
trV = tr[svmask].reset_index(drop=True)
assert len(trV) == nsv, f"truth SV {len(trV)} != panel SV {nsv} -- order mismatch"
truth_af = trV["truth_af"].values
sv_size = (trV["alt_len"] - trV["ref_len"]).abs().values

# ---- 3. kMate AF in SV-mask order -------------------------------------------
print("[3] reading kMate global ...", file=sys.stderr)
km = pd.read_csv(KMATE, sep="\t")
km = km[np.abs(km.alt_len - km.ref_len) >= SVLEN].reset_index(drop=True)
assert len(km) == nsv, f"kMate SV {len(km)} != panel SV {nsv}"
kmate_af = km["alt_freq"].values

# ---- 4. vg AF (call -v AD) placed at svidx by exact sequence -----------------
print("[4] reading vg call -v ...", file=sys.stderr)
vg_af = np.full(nsv, np.nan)
with open(VG_VCF) as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        chrom, pos, ref, alts = f[0], int(f[1]), f[3], f[4].split(",")
        fmt = f[8].split(":")
        if "AD" not in fmt:
            continue
        ad = f[9].split(":")[fmt.index("AD")].split(",")
        try:
            ad = [int(x) for x in ad]
        except ValueError:
            continue
        ref_ad = ad[0] if ad else 0
        for i, alt in enumerate(alts):
            if alt in (".", "*", "<NON_REF>") or abs(len(alt) - len(ref)) < SVLEN:
                continue
            sv = key2idx.get(seq_key(chrom, pos, ref, alt))
            if sv is None:
                continue
            alt_ad = ad[i + 1] if i + 1 < len(ad) else 0
            dp = ref_ad + alt_ad
            vg_af[sv] = (alt_ad / dp) if dp > 0 else np.nan
print(f"    vg placed {np.isfinite(vg_af).sum():,}/{nsv:,} SVs", file=sys.stderr)

# ---- 5. stratify ------------------------------------------------------------
df = pd.DataFrame(dict(svidx=np.arange(nsv), pos=pos_of_idx, mult=mult,
                       sv_size=sv_size, n_called=n_called,
                       truth=truth_af, vg=vg_af, kmate=kmate_af))

def block(name, sub):
    print(f"\n### {name}  (n={len(sub):,})")
    for tool in ("kmate", "vg"):
        mt = metrics(sub["truth"].values, sub[tool].values)
        print(f"  {tool:6s} n={mt['n']:>6} MAE={mt['MAE']:.4f} R2={mt['R2']:>7.4f} "
              f"r={mt['r']:.4f} slope={mt['slope']:.3f}")

print("=" * 64)
block("ALL SVs", df)
mult_bins = [(1, 1, "biallelic (1 ALT)"), (2, 3, "2-3 ALTs"),
             (4, 9, "4-9 ALTs"), (10, 10**9, "10+ ALTs")]
rows = []
for lo, hi, lab in mult_bins:
    sub = df[(df["mult"] >= lo) & (df["mult"] <= hi)]
    block(f"multiplicity {lab}", sub)
    for tool in ("kmate", "vg"):
        mt = metrics(sub["truth"].values, sub[tool].values)
        rows.append(dict(stratum=lab, tool=tool, **mt))

# also: vg on biallelic + adequately covered (apples to Ash regime)
sub = df[(df["mult"] == 1) & (df["n_called"] >= 0.5 * F)]
block("biallelic & panel-call-rate>=50% (Ash-like clean regime)", sub)

pd.DataFrame(rows).to_csv(str(OUTPFX) + "_metrics.tsv", sep="\t", index=False)

# ---- 6. figure --------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
def scat(a, x, y, title):
    mok = np.isfinite(x) & np.isfinite(y)
    a.hexbin(x[mok], y[mok], gridsize=60, bins="log", cmap="viridis", mincnt=1)
    a.plot([0, 1], [0, 1], "r--", lw=1)
    mt = metrics(x, y)
    a.set_title(f"{title}\nr={mt['r']:.3f}  R$^2$={mt['R2']:.3f}  slope={mt['slope']:.2f}  n={mt['n']:,}")
    a.set_xlabel("truth alt-AF"); a.set_ylabel("estimated alt-AF")
    a.set_xlim(0, 1); a.set_ylim(0, 1)

scat(ax[0], df["truth"].values, df["kmate"].values, "kMate (all SVs)")
scat(ax[1], df["truth"].values, df["vg"].values, "vg call -v (all SVs)")
# panel C: vg accuracy vs multiplicity
labs = [l for _, _, l in mult_bins]
vg_r = [metrics(df[(df["mult"] >= lo) & (df["mult"] <= hi)]["truth"].values,
                df[(df["mult"] >= lo) & (df["mult"] <= hi)]["vg"].values)["r"]
        for lo, hi, _ in mult_bins]
km_r = [metrics(df[(df["mult"] >= lo) & (df["mult"] <= hi)]["truth"].values,
                df[(df["mult"] >= lo) & (df["mult"] <= hi)]["kmate"].values)["r"]
        for lo, hi, _ in mult_bins]
xb = np.arange(len(labs))
ax[2].plot(xb, km_r, "o-", label="kMate", color="#2ca02c")
ax[2].plot(xb, vg_r, "s-", label="vg call -v", color="#7b3294")
ax[2].set_xticks(xb); ax[2].set_xticklabels(labs, rotation=20, ha="right", fontsize=8)
ax[2].set_ylabel("Pearson r vs truth"); ax[2].set_ylim(0, 1)
ax[2].set_title("accuracy vs locus allele-multiplicity")
ax[2].legend(); ax[2].grid(alpha=0.3)
fig.suptitle(f"vg-SV underperformance is driven by multiallelic cactus snarls  —  {POOL}", y=1.02)
fig.tight_layout()
fig.savefig(str(OUTPFX) + ".png", dpi=130, bbox_inches="tight")
print(f"\n[fig] {OUTPFX}.png")
print(f"[tsv] {OUTPFX}_metrics.tsv")

# readable single-panel truth-vs-est scatters (kMate, vg) ---------------------
for tag, col, title in [("kmate", "kmate", "kMate"), ("vg", "vg", "vg (construct + call -v)")]:
    f1, a1 = plt.subplots(figsize=(7, 7))
    x = df["truth"].values; y = df[col].values
    mok = np.isfinite(x) & np.isfinite(y)
    hb = a1.hexbin(x[mok], y[mok], gridsize=80, bins="log", cmap="viridis", mincnt=1)
    a1.plot([0, 1], [0, 1], "r--", lw=1.2)
    mt = metrics(x, y)
    a1.set_title(f"{title} — SV alt-AF  ({POOL.split('_')[0]}, n80 g0)\n"
                 f"r={mt['r']:.3f}  R²={mt['R2']:.3f}  MAE={mt['MAE']:.4f}  "
                 f"slope={mt['slope']:.2f}  n={mt['n']:,}", fontsize=12)
    a1.set_xlabel("truth alt-AF"); a1.set_ylabel(f"{title} estimated alt-AF")
    a1.set_xlim(0, 1); a1.set_ylim(0, 1)
    f1.colorbar(hb, ax=a1, label="log10(count)")
    f1.tight_layout()
    f1.savefig(f"{OUTPFX}__{tag}_scatter.png", dpi=140, bbox_inches="tight")
    print(f"[fig] {OUTPFX}__{tag}_scatter.png")
