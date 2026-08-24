#!/usr/bin/env python3
"""Are climate-adaptive haploblocks enriched for structural variants?

Joins the per-block SV landscape to the multisite founder-GWAS scores (JOINT /
GLOBAL / CLIMATE, per clq0.9 block) and tests SV enrichment in the top-ranked
("adaptive") blocks against a SIZE-MATCHED permutation null — because block size
(n_kept) drives both SV presence (0->47%) and GWAS power, so it MUST be matched.

Primary metric = has_sv (with tight n_kept matching). Secondary = sv_frac, n_sv.
Run in `basic` env. Deterministic (seed=0)."""
from pathlib import Path
import numpy as np, pandas as pd

OUT = Path("analysis/grenenet_gea/sv_adaptive/results")
NPERM = 10000
RNG = np.random.default_rng(0)
CONTRASTS = [("JOINT", "p_joint"), ("GLOBAL", "p_global"), ("CLIMATE", "p_clim")]
TOP = [0.005, 0.01, 0.02]          # top fractions defining "adaptive"
# non-SNP = small-indel + SV = all length-changing variation (~23% of records, far
# more abundant than SV alone). frac is size-robust; also break out indel vs sv.
METRICS = ["snp_frac", "nonsnp_frac", "has_nonsnp", "n_nonsnp", "indel_frac", "sv_frac"]
# fine n_kept bins (dense at the low end, where most blocks live)
EDGES = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 70, 100, 150, 250, 10**9]

g = pd.read_csv(OUT.parent / "hapfreq/multisite_founder_gwas_clq90_pc1.csv")
gc = g.groupby("unit").agg(p_joint=("p_joint", "min"), p_global=("p_global", "min"),
                           p_clim=("p_clim", "min")).reset_index()
L = pd.read_csv(OUT / "sv_landscape_clq0.9.csv")
df = gc.merge(L, left_on="unit", right_on="block_id", how="inner").reset_index(drop=True)
# non-SNP (indel+SV) derived metrics
df["n_nonsnp"] = df.n_indel + df.n_sv
df["nonsnp_frac"] = df.n_nonsnp / df.n_kept
df["indel_frac"] = df.n_indel / df.n_kept
df["snp_frac"] = df.n_snp / df.n_kept                  # the ~77% background, for a direct class bar
df["has_nonsnp"] = (df.n_nonsnp > 0).astype(int)
df["bin"] = pd.cut(df.n_kept, EDGES, right=False, labels=False)
# drop blocks below the size-bin floor (n_kept<2): at the common-SV floor (MAC>=12) some
# blocks retain 0-1 kept variants -> NaN bin -> unmatchable. Same handling as
# sv_frac_markermatched.py; harmless at MAC>=2 (min n_kept was already 2).
n_pre = len(df)
df = df[df.bin.notna()].reset_index(drop=True); df["bin"] = df.bin.astype(int)
if len(df) < n_pre:
    print(f"dropped {n_pre - len(df):,} blocks with n_kept<2 (below size-bin floor)")
bin_members = {b: df.index[df.bin == b].to_numpy() for b in df.bin.unique()}
print(f"{len(df):,} blocks joined | non-SNP: {int(df.has_nonsnp.sum()):,} blocks carry >=1 "
      f"({100*df.has_nonsnp.mean():.1f}%); non-SNP = {100*df.n_nonsnp.sum()/df.n_kept.sum():.1f}% of records "
      f"| SV-only: {100*df.has_sv.mean():.1f}% of blocks")

# diagnostic: is each contrast's score just a size proxy?
print("\nadaptation<->size coupling  corr(-log10 p , log10 n_kept):")
for name, pc in CONTRASTS:
    print(f"  {name:8s}: {np.corrcoef(-np.log10(df[pc].clip(1e-300)), np.log10(df.n_kept))[0,1]:+.3f}")


def matched_perm(adaptive_idx, metric):
    vals = df[metric].to_numpy()
    obs = vals[adaptive_idx].mean()
    # per adaptive block, draw NPERM controls from its own size bin
    null = np.empty((len(adaptive_idx), NPERM))
    for k, i in enumerate(adaptive_idx):
        pool = vals[bin_members[df.bin.iat[i]]]
        null[k] = RNG.choice(pool, NPERM)
    nm = null.mean(0)
    p = (1 + (nm >= obs).sum()) / (NPERM + 1)
    return obs, nm.mean(), obs / nm.mean() if nm.mean() > 0 else np.nan, p


rows = []
for name, pc in CONTRASTS:
    order = df[pc].to_numpy().argsort()          # smallest p = most adaptive
    for frac in TOP:
        n = int(round(frac * len(df)))
        adaptive = order[:n]
        med_nk = df.n_kept.iloc[adaptive].median()
        for metric in METRICS:
            obs, null, ratio, p = matched_perm(adaptive, metric)
            rows.append(dict(contrast=name, top=f"{frac:.1%}", n_adaptive=n,
                             med_nkept=med_nk, metric=metric,
                             observed=round(obs, 4), null=round(null, 4),
                             enrich=round(ratio, 2), p_perm=round(p, 4)))

res = pd.DataFrame(rows)
res.to_csv(OUT / "sv_enrichment.csv", index=False)
pd.set_option("display.width", 160, "display.max_rows", 200)
print("\n== SIZE-MATCHED SV enrichment (observed vs size-matched null) ==")
print(res.to_string(index=False))
print(f"\n-> {OUT/'sv_enrichment.csv'}")

# headline: non-SNP fraction at top 1%
print("\nHEADLINE (non-SNP fraction, top 1%, size-matched):")
for name, _ in CONTRASTS:
    r = res[(res.contrast == name) & (res.top == "1.0%") & (res.metric == "nonsnp_frac")].iloc[0]
    verdict = "ENRICHED" if (r.p_perm < 0.05 and r.enrich > 1) else ("depleted" if r.enrich < 1 and r.p_perm > 0.95 else "n.s.")
    print(f"  {name:8s}: adaptive non-SNP frac {r.observed:.3f} vs "
          f"{r.null:.3f} size-matched  (x{r.enrich}, p={r.p_perm})  {verdict}")
