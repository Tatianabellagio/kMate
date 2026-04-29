#!/usr/bin/env python3
"""Generate scatter plots for the seedmix 3-way comparison."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path("/home/tbellagio/scratch/hapfire_sv")
df = pd.read_csv(ROOT / "results" / "seedmix_3way_comparison.tsv.gz", sep="\t")
sample = "SEEDMIX_S1"
sub = df[df["sample_id"] == sample].copy()
sub_freqk = sub.dropna(subset=["freqk_af"]).copy()
print(f"S1 panel SVs: {len(sub)}; with freqk: {len(sub_freqk)}")

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Row 1: full panel (hapfire_proj vs truth_panel) — no freqk subset
ax = axes[0, 0]
ax.scatter(sub["truth_panel"], sub["hapfire_proj_af"], s=4, alpha=0.4, c="#1f77b4")
m = float(max(sub["truth_panel"].max(), sub["hapfire_proj_af"].max())) * 1.05
ax.plot([0, m], [0, m], "k--", lw=0.5)
r = sub[["truth_panel", "hapfire_proj_af"]].corr().iloc[0, 1]
mae = (sub["hapfire_proj_af"] - sub["truth_panel"]).abs().mean()
ax.set_xlabel("truth_panel (recipe @ 80 founders)")
ax.set_ylabel("hapFIRE projection")
ax.set_title(f"hapFIRE proj vs truth_panel  (n={len(sub):,}; r={r:.3f}; MAE={mae:.4f})")

# Row 1: hapfire_proj vs freqk over the freqk-overlap
ax = axes[0, 1]
ax.scatter(sub_freqk["freqk_af"], sub_freqk["hapfire_proj_af"], s=4, alpha=0.4, c="#ff7f0e")
m = max(sub_freqk["freqk_af"].max(), sub_freqk["hapfire_proj_af"].max()) * 1.05
ax.plot([0, m], [0, m], "k--", lw=0.5)
r = sub_freqk[["freqk_af", "hapfire_proj_af"]].corr().iloc[0, 1]
ax.set_xlabel("freqk AF")
ax.set_ylabel("hapFIRE projection")
ax.set_title(f"hapFIRE proj vs freqk  (n={len(sub_freqk):,}; r={r:.3f})")

# Row 1: freqk vs truth_panel
ax = axes[0, 2]
ax.scatter(sub_freqk["truth_panel"], sub_freqk["freqk_af"], s=4, alpha=0.4, c="#2ca02c")
m = max(sub_freqk["truth_panel"].max(), sub_freqk["freqk_af"].max()) * 1.05
ax.plot([0, m], [0, m], "k--", lw=0.5)
r = sub_freqk[["truth_panel", "freqk_af"]].corr().iloc[0, 1]
ax.set_xlabel("truth_panel (recipe @ 80 founders)")
ax.set_ylabel("freqk AF")
ax.set_title(f"freqk vs truth_panel  (n={len(sub_freqk):,}; r={r:.3f})")

# Row 2: same but split by var_type
for j, (vt, color) in enumerate([("DEL", "#d62728"), ("INS", "#9467bd"), ("ALL", "#7f7f7f")]):
    ax = axes[1, j]
    if vt == "ALL":
        s = sub_freqk
    else:
        s = sub_freqk[sub_freqk["var_type"] == vt]
    ax.scatter(s["truth_panel"], s["hapfire_proj_af"], s=3, alpha=0.4, c=color, label="hapFIRE proj")
    ax.scatter(s["truth_panel"], s["freqk_af"], s=3, alpha=0.2, c="grey", label="freqk")
    m = float(max(s["truth_panel"].max(), max(s["hapfire_proj_af"].max(), s["freqk_af"].max()))) * 1.05
    ax.plot([0, m], [0, m], "k--", lw=0.5)
    ax.set_xlabel("truth_panel")
    ax.set_ylabel("estimate")
    r_h = float(s[["truth_panel","hapfire_proj_af"]].corr().iloc[0,1])
    r_f = float(s[["truth_panel","freqk_af"]].corr().iloc[0,1])
    ax.set_title(f"{vt}  n={len(s):,}  r_hapf={r_h:.3f}  r_freqk={r_f:.3f}")
    ax.legend(loc="lower right", fontsize=7)

fig.suptitle(f"SEEDMIX_S1: hapFIRE projection vs freqk vs panel-restricted truth\n"
             f"truth_panel = Σ_{{e ∈ panel∩recipe}} seed_prop(e) · G_SV(e, v); panel covers ~17.7% of seedmix mass",
             fontsize=11)
fig.tight_layout()
out = ROOT / "results" / "seedmix_S1_comparison.png"
fig.savefig(out, dpi=120)
print(f"wrote: {out}")

# also a per-sample summary correlation bar chart across S1..S8
samples = sorted(df["sample_id"].unique())
rows = []
for sid in samples:
    s_all = df[df["sample_id"] == sid]
    s_frq = s_all.dropna(subset=["freqk_af"])
    rows.append({
        "sample_id": sid,
        "r_hapf_vs_truth_all": float(s_all[["hapfire_proj_af","truth_panel"]].corr().iloc[0,1]),
        "r_freqk_vs_truth": float(s_frq[["freqk_af","truth_panel"]].corr().iloc[0,1]),
        "r_hapf_vs_freqk": float(s_frq[["hapfire_proj_af","freqk_af"]].corr().iloc[0,1]),
    })
summary = pd.DataFrame(rows)
print(summary.round(3).to_string(index=False))
summary.to_csv(ROOT / "results" / "seedmix_per_sample_summary.tsv", sep="\t", index=False)

fig2, ax = plt.subplots(figsize=(9, 4))
x = np.arange(len(samples))
w = 0.27
ax.bar(x - w, summary["r_hapf_vs_truth_all"], w, label="hapFIRE-proj vs truth_panel (all SVs)", color="#1f77b4")
ax.bar(x,     summary["r_freqk_vs_truth"],   w, label="freqk vs truth_panel (overlap)",        color="#2ca02c")
ax.bar(x + w, summary["r_hapf_vs_freqk"],    w, label="hapFIRE-proj vs freqk (overlap)",       color="#ff7f0e")
ax.set_xticks(x); ax.set_xticklabels(samples, rotation=30, ha="right")
ax.set_ylim(0, 1.05); ax.set_ylabel("Pearson r")
ax.legend(loc="lower right", fontsize=8)
ax.set_title("Per-sample SV freq correlations on the 8 SEEDMIX replicates")
ax.grid(True, axis="y", alpha=0.3)
fig2.tight_layout()
out2 = ROOT / "results" / "seedmix_correlation_bars.png"
fig2.savefig(out2, dpi=120)
print(f"wrote: {out2}")
