#!/usr/bin/env python3
"""Are POOL-SEQ TEMPORAL climate-adaptive HAPLOTYPES enriched for structural variants?

GEA-realm counterpart of sv_enrichment.py (which used the founder GWAS). The adaptation
axis here is the clq0.9 Pipeline B pool-temporal climate GEA at the **HapFM-haplotype**
level (`hap_gea.csv`): every HapFM haplotype — a founder cluster within an r2>=0.9 block —
gets ONE temporal selection coefficient `s_mean` and ONE climate gradient `beta1`
(site-permutation p_two / p_one).  This is the user's "one selection coefficient per
haploblock, then test the temperature gradient" — NO within-block Sidak aggregation.

Each haplotype inherits the SV landscape of its parent clq0.9 block (sv_landscape_clq0.9).
"Adaptive" = top-X% haplotypes by climate p (or by |s_mean| for the any-direction selection
contrast).  Same SIZE-MATCHED permutation null as the GWAS pass — block `n_kept` drives BOTH
SV presence AND GEA power, so each adaptive haplotype is compared only to haplotypes whose
block sits in the same size bin.  Run in `basic`.  Deterministic (seed=0).

Contrasts:
  GEA_climate    two-sided climate association (p_two)  — the headline "climate-adaptive"
  GEA_warm       up-in-warm directional         (p_one)
  GEA_selection  strongest temporal mover, any direction (|s_mean|) — GEA analog of GWAS JOINT
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.special import erfcinv

OUT = Path("analysis/grenenet_gea/sv_adaptive/results")
GEADIR = Path("results/grenenet_gea/hapfreq_clq90/pipelineB_varlen")
NPERM = 10000
RNG = np.random.default_rng(0)
TOP = [0.005, 0.01, 0.02]
METRICS = ["sv_frac", "has_sv", "n_sv", "indel_frac", "nonsnp_frac"]
EDGES = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 70, 100, 150, 250, 10**9]

# ---- SV landscape per clq0.9 block -------------------------------------------------
L = pd.read_csv(OUT / "sv_landscape_clq0.9.csv")
L["nonsnp_frac"] = (L.n_indel + L.n_sv) / L.n_kept
L["indel_frac"] = L.n_indel / L.n_kept

# ---- per-HapFM-haplotype climate GEA ----------------------------------------------
hg = pd.read_csv(GEADIR / "hap_gea.csv")
hg["unit"] = hg.chrom + ":" + hg.unit_start.astype(str) + "-" + hg.unit_end.astype(str)
# testable set == the GEA test's own: covered + panel_freq in [0.05,0.95], and keep the
# k-1 INDEPENDENT (non-reference) haplotypes per block (drop the single highest-panel_freq
# member = the 'rest' reference; k=2 -> keep the minor only).  Matches build_hap_gea.
hg = hg[hg.covered & hg.panel_freq.between(0.05, 0.95)].copy()
# drop the single highest-panel_freq member per unit (the 'rest' reference), keeping
# units with only one member intact. (pandas 3.0: groupby(...).apply drops the grouping
# column, so build the drop-index explicitly instead of returning trimmed groups.)
_sizes = hg.groupby("unit").size()
_drop = hg.groupby("unit")["panel_freq"].idxmax()
hg = hg.drop(index=_drop[_sizes > 1].to_numpy())

df = hg.merge(L, left_on="unit", right_on="block_id", how="inner").reset_index(drop=True)
df["bin"] = pd.cut(df.n_kept, EDGES, right=False, labels=False)
# drop haplotypes whose block has n_kept<2 (NaN bin -> unmatchable). At the common-SV floor
# (MAC>=12) some blocks retain 0-1 kept variants; harmless at MAC>=2 (min n_kept was 2).
n_pre = len(df)
df = df[df.bin.notna()].reset_index(drop=True); df["bin"] = df.bin.astype(int)
if len(df) < n_pre:
    print(f"dropped {n_pre - len(df):,} haplotypes in blocks with n_kept<2 (below size-bin floor)")
bin_members = {b: df.index[df.bin == b].to_numpy() for b in df.bin.unique()}

# Genomic inflation + permutation-floor tie mass: with 5000 perms the p_two floor is
# 1/5001, and a large fraction of haplotypes tie there (their observed |beta1| beats ALL
# permutations). So p_two alone cannot rank the top few % — it must be broken by the
# CONTINUOUS climate effect size |beta1|.  Absolute significance is not claimed (the
# N~19-site pseudoreplication wall); this is a RELATIVE, ranking-based contrast.
nfloor = (df.p_two == df.p_two.min()).sum()
lam = np.median((np.sqrt(2) * erfcinv(df.p_two.to_numpy()))**2) / 0.4549    # genomic inflation
print(f"{len(df):,} testable HapFM haplotypes carry a clq90 climate-GEA score, "
      f"spanning {df.unit.nunique():,} clq0.9 blocks "
      f"({100*df.unit.nunique()/len(L):.0f}% of the {len(L):,} landscape blocks)")
print(f"  haplotypes whose block has >=1 SV: {int(df.has_sv.sum()):,} ({100*df.has_sv.mean():.1f}%)")
print(f"  permutation floor p={df.p_two.min():.5f}: {nfloor:,} haplotypes tied ({100*nfloor/len(df):.2f}%) "
      f"| genomic inflation lambda={lam:.2f}  -> rank the floor by |beta1|, no absolute p claimed")

# ordering: adaptive-first.  climate = (p_two asc, |beta1| desc) so the floor-tied set is
# ordered by climate-gradient magnitude; warm = (p_one asc, beta1 desc, most up-in-warm);
# selection = |s_mean| desc (continuous temporal-mover magnitude, GEA analog of GWAS JOINT).
ORDER = {
    "GEA_climate":   np.lexsort((-df.beta1.abs().to_numpy(), df.p_two.to_numpy())),
    "GEA_warm":      np.lexsort((-df.beta1.to_numpy(),        df.p_one.to_numpy())),
    "GEA_selection": np.argsort(-df.s_mean.abs().to_numpy(), kind="stable"),
}
CONTRASTS = list(ORDER)
print("\nadaptation<->size coupling  corr(effect-size, log10 n_kept)  [why size-matching is load-bearing]:")
for nm, col in [("GEA_climate", df.beta1.abs()), ("GEA_warm", df.beta1.abs()),
                ("GEA_selection", df.s_mean.abs())]:
    print(f"  {nm:13s}: corr(effect-size, log10 n_kept) = {np.corrcoef(col, np.log10(df.n_kept))[0,1]:+.3f}")


def matched_perm(adaptive_idx, metric):
    vals = df[metric].to_numpy()
    obs = vals[adaptive_idx].mean()
    null = np.empty((len(adaptive_idx), NPERM))
    for k, i in enumerate(adaptive_idx):
        null[k] = RNG.choice(vals[bin_members[df.bin.iat[i]]], NPERM)
    nm = null.mean(0)
    p = (1 + (nm >= obs).sum()) / (NPERM + 1)
    return obs, nm.mean(), (obs / nm.mean() if nm.mean() > 0 else np.nan), p


rows = []
for name in CONTRASTS:
    order = ORDER[name]                             # adaptive-first (see ORDER above)
    for frac in TOP:
        n = int(round(frac * len(df)))
        adaptive = order[:n]
        nblk = df.unit.iloc[adaptive].nunique()
        med_nk = df.n_kept.iloc[adaptive].median()
        for metric in METRICS:
            obs, null, ratio, p = matched_perm(adaptive, metric)
            rows.append(dict(contrast=name, top=f"{frac:.1%}", n_adaptive=n, n_blocks=nblk,
                             med_nkept=med_nk, metric=metric,
                             observed=round(obs, 4), null=round(null, 4),
                             enrich=round(ratio, 2), p_perm=round(p, 4)))

res = pd.DataFrame(rows)
res.to_csv(OUT / "sv_enrichment_gea.csv", index=False)
pd.set_option("display.width", 180, "display.max_rows", 200)
print("\n== SIZE-MATCHED SV enrichment vs POOL-TEMPORAL climate GEA (clq0.9 HapFM haplotypes) ==")
print(res.to_string(index=False))
print(f"\n-> {OUT/'sv_enrichment_gea.csv'}")

print("\nHEADLINE (SV fraction, top 1%, size-matched):")
for name in CONTRASTS:
    r = res[(res.contrast == name) & (res.top == "1.0%") & (res.metric == "sv_frac")].iloc[0]
    verdict = ("ENRICHED" if (r.p_perm < 0.05 and r.enrich > 1)
               else "depleted" if (r.enrich < 1 and r.p_perm > 0.95) else "n.s.")
    print(f"  {name:13s}: adaptive SV frac {r.observed:.4f} vs {r.null:.4f} size-matched "
          f"(x{r.enrich}, p={r.p_perm}, {int(r.n_adaptive)} haps / {int(r.n_blocks)} blocks)  {verdict}")
